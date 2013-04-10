#include "InterHistThrustFunctor.hh"

__constant__ fptype* dev_base_interhists[100]; // Multiple histograms for the case of multiple PDFs
#define OBS_CODE 4242424242
// This number is presumably so high that it will never collide
// with an actual parameter index. It indicates that this dimension
// is an event observable.

// dev_powi is implemented in SmoothHistogramThrustFunctor.cu. 

__device__ fptype device_InterHistogram (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure is
  // nP totalHistograms (idx1 limit1 step1 bins1) (idx2 limit2 step2 bins2) nO o1 o2
  // where limit and step are indices into functorConstants. 

  int numVars = (indices[0] - 1) / 4; 
  int globalBin = 0; 
  int previous = 1; 
  int myHistogramIndex = indices[1]; 
  fptype binDistances[10]; // Ten dimensions should be more than enough! 
  // Distance from bin center in units of bin width in each dimension. 

  unsigned int observablesSeen = 0; 
  for (int i = 0; i < numVars; ++i) { 
    fptype currVariable = 0; 
    unsigned int varIndex = indices[2 + 4*i];
    if (varIndex == OBS_CODE) {
      // Interpret this number as observable index. 
      // Notice that this if does not cause a fork 
      // - all threads will hit the same index and 
      // make the same decision. 
      currVariable = evt[indices[indices[0] + 2 + observablesSeen++]]; 
    }
    else {
      // Interpret as parameter index. 
      currVariable = p[varIndex]; 
    }

    int lowerBoundIdx   = 3 + 4*i; 
    fptype lowerBound   = functorConstants[indices[lowerBoundIdx + 0]];
    fptype step         = functorConstants[indices[lowerBoundIdx + 1]];

    currVariable   -= lowerBound;
    currVariable   /= step; 

    int localBin    = (int) FLOOR(currVariable); 
    binDistances[i] = currVariable - localBin - fptype(0.5); 
    globalBin      += previous * localBin; 
    previous       *= indices[lowerBoundIdx + 2];

    if (0 == threadIdx.x + blockIdx.x)
      printf("Variable %i: %f %f %i\n", i, currVariable, currVariable*step + lowerBound, localBin);
  }

  fptype* myHistogram = dev_base_interhists[myHistogramIndex];
  fptype ret = 0; 

  //------------------ 
  //     |     |     |   
  //  3  |  4  |  5  |    
  //     |     |     |
  //------------------
  //    x|     |     |
  //  0  |  1  |  2  |
  //     |     |     |
  //------------------

  fptype totalWeight = 0; 
  int totalBins = dev_powi(3, numVars); 
  for (int i = 0; i < totalBins; ++i) {
    int currBin = globalBin; 
    int localPrevious = 1; 
    int trackingBin = globalBin; 
    bool offSomeAxis = false; 
    fptype currentWeight = 0; 
    // Loop over vars to get offset for each one. 
    for (int v = 0; v < numVars; ++v) {
      int localNumBins = indices[4*(v+1) + 1];
      int offset = ((i / dev_powi(3, v)) % 3) - 1; 
      
      currBin += offset * localPrevious; 
      localPrevious *= localNumBins; 
      
      int currVarBin = trackingBin % localNumBins; 
      trackingBin /= localNumBins; 
      if (currVarBin + offset < 0) offSomeAxis = true;
      if (currVarBin + offset >= localNumBins) offSomeAxis = true;

      fptype currDist = binDistances[v];
      currDist -= offset; 
      currentWeight += currDist*currDist;
      if (0 == threadIdx.x + blockIdx.x)
	printf("%i, %i: %f %f %f %i %s\n", i, v, currDist, binDistances[v], currentWeight, offset, offSomeAxis ? "off" : "on"); 
    }

    // Only interpolate the four closest boxes (in two dimensions; more in three dimensions). 
    currentWeight = currentWeight > 0 ? (currentWeight <= SQRT((fptype) numVars) ? 1 / SQRT(currentWeight) : 0) : 0; 
    fptype currentEntry = offSomeAxis ? 0 : myHistogram[currBin];
    ret += currentWeight * currentEntry;
    totalWeight += currentWeight;

    if (0 == threadIdx.x + blockIdx.x) 
      printf("Adding bin content %i %f with weight %f for total %f.\n", currBin, currentEntry, currentWeight, ret);
  }

  if (0 == threadIdx.x + blockIdx.x)
    printf("%f %f %f %i %f\n", ret, totalWeight, evt[0], indices[6], p[indices[6]]);

  ret /= totalWeight;
  return ret; 
}

__device__ device_function_ptr ptr_to_InterHistogram = device_InterHistogram; 

__host__ InterHistThrustFunctor::InterHistThrustFunctor (std::string n, 
							 BinnedDataSet* x, 
							 std::vector<Variable*> params, 
							 std::vector<Variable*> obses)
  : ThrustPdfFunctor(0, n) 
  , numVars(x->numVariables())
{
  int numConstants = 2*numVars;
  registerConstants(numConstants);
  static unsigned int totalHistograms = 0; 
  host_constants = new fptype[numConstants]; 
  totalEvents = 0; 

  std::vector<unsigned int> pindices;
  pindices.push_back(totalHistograms); 

  int varIndex = 0; 
  for (varConstIt var = x->varsBegin(); var != x->varsEnd(); ++var) {
    if (std::find(obses.begin(), obses.end(), *var) != obses.end()) {
      registerObservable(*var);
      pindices.push_back(OBS_CODE); 
    }
    else {
      pindices.push_back(registerParameter(*var));
    }

    pindices.push_back(cIndex + 2*varIndex + 0);
    pindices.push_back(cIndex + 2*varIndex + 1);
    pindices.push_back((*var)->numbins);

    // NB, do not put cIndex here, it is accounted for by the offset in cudaMemcpyToSymbol below. 
    host_constants[2*varIndex + 0] = (*var)->lowerlimit; 
    host_constants[2*varIndex + 1] = ((*var)->upperlimit - (*var)->lowerlimit) / (*var)->numbins; 
    varIndex++; 
  }

  unsigned int numbins = x->getNumBins(); 
  thrust::host_vector<fptype> host_histogram; 
  for (unsigned int i = 0; i < numbins; ++i) {
    fptype curr = x->getBinContent(i);
    host_histogram.push_back(curr);
    totalEvents += curr; 
  }
  cudaMemcpyToSymbol(functorConstants, host_constants, numConstants*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice); 

  dev_base_histogram = new thrust::device_vector<fptype>(host_histogram);  
  static fptype* dev_address[1];
  dev_address[0] = (&((*dev_base_histogram)[0])).get();
  cudaMemcpyToSymbol(dev_base_interhists, dev_address, sizeof(fptype*), totalHistograms*sizeof(fptype*), cudaMemcpyHostToDevice); 
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_InterHistogram, sizeof(void*));
  initialise(pindices); 

  totalHistograms++; 
}
/*
__host__ fptype InterHistThrustFunctor::normalise () const {
  fptype ret = thrust::reduce(dev_base_histogram->begin(), dev_base_histogram->end()); 

  for (unsigned int varIndex = 0; varIndex < numVars; ++varIndex) {
    ret *= host_constants[2*varIndex + 1]; // Bin size cached by constructor. 
  }

  ret = 50*host_params[0]; 
  host_normalisation[parameters] = 1.0/ret; 

  std::cout << getName() << " found integral " << ret << " with parameter " << host_params[0] << std::endl; 
  return ret; 
}
*/


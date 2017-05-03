#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

FCN::FCN(Params& params) : params_(&params) {
    host_callnumber = 0;
    
    // Verify that all varaibles need to be recached
    for(Variable* var : params_->vars_)
        var->unchanged_ = false;
    
}

double FCN::operator()(const std::vector<double>& pars) const {
    
    // Translate from Minuit indexing to GooFit indexing
    std::vector<double> gooPars;
    gooPars.resize(max_index(params_->vars_)+1);
    
    for(Variable* var : params_->vars_) {
        var->unchanged_ = var->value == pars.at(var->getFitterIndex());
        gooPars.at(var->getIndex()) = pars.at(var->getFitterIndex());
    }

    params_->pdf_->copyParams(gooPars);
    double nll = params_->pdf_->calculateNLL();
    host_callnumber++;

    return nll;
}

// Get the number of variable parameters
Params* FCN::GetParams() {
    return params_;
}
    
}
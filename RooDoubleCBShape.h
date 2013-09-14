#ifndef ROO_CB_SHAPE
#define ROO_CB_SHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooDoubleCBShape : public RooAbsPdf {
public:
  RooDoubleCBShape() {} ;
  RooDoubleCBShape(const char *name, const char *title, RooAbsReal& _m,
	     RooAbsReal& _m0, RooAbsReal& _sigma,
	     RooAbsReal& _alphaL, RooAbsReal& _nL, // left
	     RooAbsReal& _alphaR, RooAbsReal& _nR);  // right

  RooDoubleCBShape(const RooDoubleCBShape& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooDoubleCBShape(*this,newname); }

  inline virtual ~RooDoubleCBShape() { }

  virtual Int_t getAnalyticalIntegral( RooArgSet& allVars,  RooArgSet& analVars, const char* rangeName=0 ) const;
  virtual Double_t analyticalIntegral( Int_t code, const char* rangeName=0 ) const;

  // Optimized accept/reject generator support
  virtual Int_t getMaxVal(const RooArgSet& vars) const ;
  virtual Double_t maxVal(Int_t code) const ;

protected:

  Double_t ApproxErf(Double_t arg) const ;

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaL;
  RooRealProxy nL;
  RooRealProxy alphaR;
  RooRealProxy nR;

  Double_t evaluate() const;

private:

  ClassDef(RooDoubleCBShape,1) // Crystal Ball lineshape PDF
};

#endif

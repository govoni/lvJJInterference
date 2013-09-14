
//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// P.d.f implementing the Crystall Ball line shape
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooDoubleCBShape.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "TMath.h"

using namespace std;

ClassImp(RooDoubleCBShape)
;

  

//_____________________________________________________________________________
Double_t RooDoubleCBShape::ApproxErf(Double_t arg) const 
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return RooMath::erf(arg);
}


static const char rcsid[] =
"$Id: RooDoubleCBShape.cxx 44982 2012-07-10 08:36:13Z moneta $";


//_____________________________________________________________________________
RooDoubleCBShape::RooDoubleCBShape(const char *name, const char *title,
		       RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigma,
	     RooAbsReal& _alphaL, RooAbsReal& _nL, // left
	     RooAbsReal& _alphaR, RooAbsReal& _nR): // right
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaL("alphaL", "AlphaL", this, _alphaL),
  nL("nL", "OrderL", this, _nL),
  alphaR("alphaR", "AlphaR", this, _alphaR),
  nR("nR", "OrderR", this, _nR)
{
}


//_____________________________________________________________________________
RooDoubleCBShape::RooDoubleCBShape(const RooDoubleCBShape& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), 
  alphaL("alphaL", this, other.alphaL), nL("nL", this, other.nL),
  alphaR("alphaR", this, other.alphaR), nR("nR", this, other.nR)
{
}


//_____________________________________________________________________________
Double_t RooDoubleCBShape::evaluate() const {

  double localSigma  = fabs (sigma);

  if( (m-m0)/localSigma > fabs(alphaR) ) // right tail
  {
    double A = pow(nR/fabs(alphaR), nR) * exp(-0.5 * alphaR*alphaR);
    double B = nR/fabs(alphaR) - fabs(alphaR);
    
    return A * pow(B + (m-m0)/localSigma, -1.*nR);
  }
  
  else if( (m-m0)/localSigma < -1.*fabs(alphaL) ) // left tail
  {
    double A = pow(nL/fabs(alphaL), nL) * exp(-0.5 * alphaL*alphaL);
    double B = nL/fabs(alphaL) - fabs(alphaL);
    
    return A * pow(B - (m-m0)/localSigma, -1.*nL);
  }
  
  else // core
  {
    return exp(-1. * (m-m0)*(m-m0) / (2*localSigma*localSigma) );
  }
  
}


//_____________________________________________________________________________
Int_t RooDoubleCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1 ;
  
  return 0;
}



//_____________________________________________________________________________
Double_t RooDoubleCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{
  return 1. ; // FIXME

//  static const double sqrtPiOver2 = 1.2533141373;
//  static const double sqrt2 = 1.4142135624;
//
//  assert(code==1);
//  double result = 0.0;
//  bool useLog = false;
//  
//  if( fabs(n-1.0) < 1.0e-05 )
//    useLog = true;
//  
//  double sig = fabs((Double_t)sigma);
//  
//  double tmin = (m.min(rangeName)-m0)/sig;
//  double tmax = (m.max(rangeName)-m0)/sig;
//  
//  if(alpha < 0) {
//    double tmp = tmin;
//    tmin = -tmax;
//    tmax = -tmp;
//  }
//
//  double absAlpha = fabs((Double_t)alpha);
//  
//  if( tmin >= -absAlpha ) {
//    result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
//                                - ApproxErf(tmin/sqrt2) );
//  }
//  else if( tmax <= -absAlpha ) {
//    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
//    double b = n/absAlpha - absAlpha;
//    
//    if(useLog) {
//      result += a*sig*( log(b-tmin) - log(b-tmax) );
//    }
//    else {
//      result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
//                                - 1.0/(TMath::Power(b-tmax,n-1.0)) );
//    }
//  }
//  else {
//    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
//    double b = n/absAlpha - absAlpha;
//    
//    double term1 = 0.0;
//    if(useLog) {
//      term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
//    }
//    else {
//      term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
//                              - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
//    }
//    
//    double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
//                                     - ApproxErf(-absAlpha/sqrt2) );
//    
//    
//    result += term1 + term2;
//  }
//  
//  return result;
}



//_____________________________________________________________________________
Int_t RooDoubleCBShape::getMaxVal(const RooArgSet& vars) const 
{
  // Advertise that we know the maximum of self for given (m0,alpha,n,sigma)
  RooArgSet dummy ;

  if (matchArgs(vars,dummy,m)) {
    return 1 ;  
  }
  return 0 ;  
}



//_____________________________________________________________________________
Double_t RooDoubleCBShape::maxVal(Int_t code) const
{
  assert(code==1) ;

  // The maximum value for given (m0,alpha,n,sigma)
  return 1.0 ;
}



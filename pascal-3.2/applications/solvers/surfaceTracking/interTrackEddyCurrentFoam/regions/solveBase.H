
#   include "ABlockEqn.H"

    jRe ==  control.omega() * sigma * AIm;
    jIm == -control.omega() * sigma * ARe;

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);

    FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

    pB == 0.5 * physicalConstant::rMu0
        * 0.5 * ( (BRe & BRe) + (BIm & BIm) );

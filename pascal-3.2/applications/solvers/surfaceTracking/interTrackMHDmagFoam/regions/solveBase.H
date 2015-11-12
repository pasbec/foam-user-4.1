
#   include "ABlockEqn.H"

    jRe ==  alpha * AIm;
    jIm == -alpha * ARe;

    BRe == fvc::curl(ARe);
    BIm == fvc::curl(AIm);

    FL == 0.5 * ( (jRe ^ BRe) + (jIm ^ BIm) );

    pB == 0.5 * rMu0
        * 0.5 * ( (BRe & BRe) + (BIm & BIm) );
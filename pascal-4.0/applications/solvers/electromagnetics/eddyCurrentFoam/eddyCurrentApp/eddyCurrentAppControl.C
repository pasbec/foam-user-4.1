/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "eddyCurrentAppControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::eddyCurrentApp::Control, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Control::read()
{
    solutionControl::read(false);

    // Read solution controls
    const dictionary& eddyCurrentDict = dict();
    nCorrEDDYCURRENT_ =
        eddyCurrentDict.lookupOrDefault<label>("nCorrectors", 1);
    nSubCorrEDDYCURRENT_ =
        eddyCurrentDict.lookupOrDefault<label>("nSubCorrectors", 2147483647);
}


bool Foam::eddyCurrentApp::Control::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    forAll (mesh_, regionI)
    {
        const dictionary& solverDict =
            mesh_[regionI].solutionDict().solverPerformanceDict();

        forAllConstIter(dictionary, solverDict, iter)
        {
            const word& variableName = iter().keyword();
            const label fieldI = applyToField(variableName);
            if (fieldI != -1)
            {
                const List<solverPerformanceData> spd = iter().stream();

                const scalar residual = spd.last().maxInitialResidual();

                checked = true;

                bool absCheck = residual < residualControl_[fieldI].absTol;
                achieved = achieved && absCheck;

                if (debug)
                {
                    Info<< algorithmName_ << " abs statistics:" << endl;

                    Info<< "    " << variableName
                        << ": res = " << residual
                        << " (" << residualControl_[fieldI].absTol << ")"
                        << endl;
                }
            }
        }
    }

    return checked && achieved;
};


void Foam::eddyCurrentApp::Control::storePrevIterFields() const
{}


bool Foam::eddyCurrentApp::Control::subCriteriaSatisfied()
{
    // no checks on first sub-iteration - nothing has been calculated yet
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    forAll (mesh_, regionI)
    {
        const dictionary& solverDict =
            mesh_[regionI].solutionDict().solverPerformanceDict();

        forAllConstIter (dictionary, solverDict, iter)
        {
            const word& variableName = iter().keyword();
            const label fieldI = applyToField(variableName);
            if (fieldI != -1)
            {
                const List<solverPerformanceData> spd = iter().stream();

                scalar oldOldResidual = VGREAT;
                if (spd.size() > 2)
                {
                    oldOldResidual = spd[spd.size()-3].maxInitialResidual();
                }

                scalar oldResidual = GREAT;
                if (spd.size() > 1)
                {
                    oldResidual = spd[spd.size()-2].maxInitialResidual();
                }
                const scalar residual = spd.last().maxInitialResidual();

                scalar oldRelative =
                    mag(oldResidual/(oldOldResidual+VSMALL) - 1.0);

                scalar relative =
                    mag(residual/(oldResidual+VSMALL) - 1.0);

                checked = true;

                scalar subRelTol =
                    subScale_ * residualControl_[fieldI].relTol;

                const bool relCheck =
                    (max(oldRelative, relative) <= subRelTol);

                achieved = achieved && relCheck;

                if (debug)
                {
                    Info<< algorithmName_ << " rel statistics:" << endl;

                    Info<< "    " << variableName
                        << ": old rel res = " << oldRelative
                        << ", rel res = " << relative
                        << " (" << subRelTol << "/"
                        << residualControl_[fieldI].relTol << ")"
                        << endl;
                }
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Control::Control
(
    regionFvMesh& mesh,
    const word& dictName
)
:
    solutionControl(mesh[Region::DEFAULT], dictName),
    mesh_(mesh),
    meshIs3D_((mesh.nGeometricD() == 3)),
    solutionDir_(meshIs3D_ ? mesh.geometricD() : -mesh.geometricD()),
    interfacePatchName_
    (
        dict().lookup("interface")
    ),
    interfacePatchLabel_
    (
        mesh_[Region::CONDUCTOR].
            boundaryMesh().findPatchID(interfacePatchName_)
    ),
    nCorrEDDYCURRENT_(-1),
    nSubCorrEDDYCURRENT_(-1),
    subCorr_(0),
    subScale_(1.0)
{
    if (interfacePatchLabel_ == -1)
    {
        FatalError
            << "Conductor patch name '"
            << interfaceName()
            << "' does not exist."
            << abort(FatalError);
    }

    read();

    Info<< nl;

    if (residualControl_.empty())
    {
        Info<< algorithmName_ << ": no residual control data found. "
            << "Calculations will employ " << nCorrEDDYCURRENT_
            << " corrector loops"
            << " and " << nSubCorrEDDYCURRENT_
            << " sub-corrector loops" << nl << endl;
    }
    else
    {
        Info<< algorithmName_ << ": max iterations = "
            << nCorrEDDYCURRENT_
            << endl
            << algorithmName_ << ": max sub-iterations = "
            << nSubCorrEDDYCURRENT_
            << endl;
        forAll(residualControl_, i)
        {
            Info<< "    field " << residualControl_[i].name << token::TAB
                << ": rel tol " << residualControl_[i].relTol
                << ", tol " << residualControl_[i].absTol
                << nl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::eddyCurrentApp::Control::loop()
{
    read();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    if (corr_ == nCorrEDDYCURRENT_ + 1)
    {
        Info<< algorithmName_ << ": not converged within "
            << nCorrEDDYCURRENT_ << " iterations" << endl;

        corr_ = 0;
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }

    if ((corr_ > 1) && criteriaSatisfied())
    {
        Info<< algorithmName_ << ": converged in " << corr_ - 1
            << " iterations" << endl;

        corr_ = 0;
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }
    else
    {
        Info<< nl;
        Info<< algorithmName_ << ": iteration " << corr_ << endl;
        storePrevIterFields();
    }

    return true;
};


Foam::dictionary Foam::eddyCurrentApp::Control::subDict
(
    label regionI,
    const word& name
) const
{
    dictionary dict = mesh_[regionI].
        solutionDict().subDict("solvers").subDict(name);

    scalar relTol = readScalar(dict.lookup("relTol"));

    scalar relTolScale = VGREAT;
    scalar progressLeft;

    if ((corr_ == 1) && (subCorr_ == 1))
    {
        progressLeft = 1.0;
    }
    else
    {
        progressLeft = 0.0;

        forAll (mesh_, regioni)
        {
            const dictionary& solverDict =
                mesh_[regionI].solutionDict().solverPerformanceDict();

            forAllConstIter (dictionary, solverDict, iter)
            {
                const word& variableName = iter().keyword();
                const label fieldI = applyToField(variableName);
                if (fieldI != -1)
                {
                    const List<solverPerformanceData> spd = iter().stream();

                    const scalar residual = spd.last().maxInitialResidual();

                    // Residual difference from target residual
                    scalar residualDiff =
                        residual-residualControl_[fieldI].absTol;

                    // Relative tolerance scaling depending on absolute tol
                    relTolScale = pow(residualControl_[fieldI].absTol, -1);

                    // Linear convergence progress left
                    progressLeft =
                        max
                        (
                            progressLeft,
                            min(max(residualDiff, 0.0), 1.0)
                        );
                }
            }
        }
    }

    // Relative tolerance decrease function (atan, scaled)
    scalar decRelTol =
        2.0 * atan(subScale_ * relTolScale * progressLeft) / M_PI;

    // Relative tolerance of sub-iteration
    scalar subRelTol = decRelTol * relTol;

    // Overwrite relTol in dictionary
    dict.set<scalar>
    (
        "relTol",
        subRelTol
    );

    if (debug)
    {
        Info<< algorithmName_ << " sub-loop:" << endl;

        Info<< "    " << name
            << " " << algorithmName_ << " iter dict " << subCorr_
            << ": progress left = " << progressLeft
            << ", rel tol scale = " << relTolScale
            << ", rel tol = " << subRelTol << "/" << relTol
            << endl;
    }

    return dict;
};


bool Foam::eddyCurrentApp::Control::subLoop()
{
    subCorr_++;

    if (debug)
    {
        Info<< algorithmName_ << " sub-loop: corr = " << subCorr_ << endl;
    }

    if ((subCorr_ > 1) && criteriaSatisfied())
    {
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }

    if (subCorr_ == nSubCorrEDDYCURRENT_ + 1)
    {
        subCorr_ = 0;
        subScale_ /= 10.0;

        Info<< algorithmName_ << " sub-loop: scale = " << subScale_ << endl;

        return false;
    }

    if ((subCorr_ > 1) && subCriteriaSatisfied() && meshIs3D())
    {
        subCorr_ = 0;
        subScale_ /= 10.0;

        Info<< algorithmName_ << " sub-loop: scale = " << subScale_ << endl;

        return false;
    }
    else
    {
        Info<< algorithmName_ << ": sub-iteration " << subCorr_ << endl;
        storePrevIterFields();
    }

    return true;
};


// ************************************************************************* //

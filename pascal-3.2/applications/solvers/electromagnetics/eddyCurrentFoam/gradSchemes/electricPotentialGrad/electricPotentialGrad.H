/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    Foam::fv::electricPotentialGrad

Description
    TODO

    NOTE: Currently it is NOT SUPPORTED or even tested to put the
          conducting region on any type of boundary

SourceFiles
    electricPotentialGrad.C

\*---------------------------------------------------------------------------*/

#ifndef electricPotentialGrad_H
#define electricPotentialGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"
#include "faceSet.H"
#include "electricPotentialGradBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class electricPotentialGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class electricPotentialGrad
:
    public fv::electricPotentialGradBase,
    public fv::gradScheme<Type>
{

private:

    // Private data

        //- Interpolation scheme
        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;

        //- Magnetic properties dictionary
        const IOdictionary& magneticProperties_;

        //- Name of the faceSet for non-conducting/conducting interface
        word nameInterface_;

        //- Name of the conductivity
        word nameConductivity_;

        //- Name of the electric scalar potential
        word nameElectricPotential_;

        //- Name of the magnetic vector potential
        word nameMagneticPotential_;


    // Private Member Functions

        //- Read dictionary data
        void readMagneticProperties();

        //- Disallow default bitwise copy construct
        electricPotentialGrad(const electricPotentialGrad&);

        //- Disallow default bitwise assignment
        void operator=(const electricPotentialGrad&);


public:

    //- Runtime type information
    TypeName("electricPotential");


    // Constructors

        //- Construct from mesh
        electricPotentialGrad
        (
            const fvMesh& mesh
        )
        :
            electricPotentialGradBase(GAUSS),
            gradScheme<Type>(mesh),
            tinterpScheme_
            (
                new linear<Type>
                (
                    mesh
                )
            ),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            )
        {
            checkBaseScheme();

            readMagneticProperties();
        }

        //- Construct from mesh and Istream
        electricPotentialGrad
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            electricPotentialGradBase(is),
            gradScheme<Type>(mesh),
            tinterpScheme_(NULL),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            )
        {
            checkBaseScheme();

            if (baseScheme() == GAUSS)
            {
                tinterpScheme_ =
                tmp<surfaceInterpolationScheme<Type> >
                (
                    surfaceInterpolationScheme<Type>::New
                    (
                        mesh,
                        is
                    )
                );
            }

            readMagneticProperties();
        }


    // Member Functions

        //- Return the gradient of the given field calculated
        //  using Gauss' theorem on the interpolated field
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;

        //- Return the BlockLduSystem corresponding to the implicit grad
        // discretization. For block coupled systems.
        virtual tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarGaussGrad.H"

#ifdef NoRepository
#   include "electricPotentialGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

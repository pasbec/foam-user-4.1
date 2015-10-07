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
#include "fluxConservativeValue.H"
#include "linear.H"

#include "faceSet.H"

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
    public fv::gradScheme<Type>
{
    // Private data

        //- Magnetic properties dictionary
        const IOdictionary& magneticProperties_;

        //- Flux conservative face value interpolation scheme
        tmp<fluxConservativeValue<Type> > tinterpScheme_;

        //- Name of the faceSet for non-conducting/conducting interface
        word nameInterface_;

        //- Name of the conductivity
        word nameConductivity_;

        //- Name of the electric scalar potential
        word nameElectricPotential_;

        //- Name of the magnetic vector potential
        word nameMagneticPotential_;


    // Private Member Functions

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
            const fvMesh& mesh,
            word nameMagneticProperties
        )
        :
            gradScheme<Type>(mesh),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            tinterpScheme_
            (
                new fluxConservativeValue<Type>
                (
                    word(magneticProperties_.lookup("conductivity"))
                )
            )
        {
            dictionary namesDict(magneticProperties_.subDict("names"));

            nameInterface_ = word(namesDict.lookup("interface"));
            nameConductivity_ = word(namesDict.lookup("conductivity"));
            nameElectricPotential_ = word(namesDict.lookup("electricPotential"));
            nameMagneticPotential_ = word(namesDict.lookup("magneticPotential"));
        }

        //- Construct from mesh and Istream
        electricPotentialGrad
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            magneticProperties_
            (
                mesh.lookupObject<IOdictionary>
                (
                    "magneticProperties"
                )
            ),
            tinterpScheme_
            (
                new fluxConservativeValue<Type>
                (
                    mesh.lookupObject<volScalarField>
                    (
                        word
                        (
                            magneticProperties_.subDict("names")
                            .lookup("conductivity")
                        )
                    )
                )
            )
        {
            dictionary namesDict(magneticProperties_.subDict("names"));

            nameInterface_ = word(namesDict.lookup("interface"));
            nameConductivity_ = word(namesDict.lookup("conductivity"));
            nameElectricPotential_ = word(namesDict.lookup("electricPotential"));
            nameMagneticPotential_ = word(namesDict.lookup("magneticPotential"));
        }


    // Member Functions

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > gradf
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
            const word& name
        ) const;

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
        tmp
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

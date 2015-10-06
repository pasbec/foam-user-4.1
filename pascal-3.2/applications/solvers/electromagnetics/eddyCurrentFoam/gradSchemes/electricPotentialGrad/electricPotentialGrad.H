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

        //- Base gradient scheme
        word baseGradSchemeName_;

        //- Interpolation scheme
        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;

        //- Name magnetic properties dictionary
        word nameMagneticProperties_;

        word nameInterface_;

        word nameScalarPotential_;

        word nameVectorPotential_;

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
            tinterpScheme_(new linear<Type>(mesh)),
            nameMagneticProperties_("magneticProperties")
        {}

        //- Construct from mesh and Istream
        electricPotentialGrad
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(NULL)
        {
            nameMagneticProperties_ = word("magneticProperties");

            if (is.eof())
            {
                FatalIOErrorIn
                (
                    "electricPotentialGrad<Type>::electricPotentialGrad(const fvMesh&, Istream&)",
                    is
                )   << "Name of magnetic vector potential not specified"
                    << exit(FatalIOError);
            }

            dictionary namesDict
            (
                mesh.lookupObject<IOdictionary>
                (
                    nameMagneticProperties_
                ).subDict("names")
            );

            nameInterface_ = word(namesDict.lookup("interface"));
            nameScalarPotential_ = word(namesDict.lookup("scalarPotential"));
            nameVectorPotential_ = word(namesDict.lookup("vectorPotential"));

            baseGradSchemeName_ = word(is);

            if (baseGradSchemeName_ != "Gauss")
            {
                FatalIOErrorIn
                (
                    "electricPotentialGrad<Type>::electricPotentialGrad(const fvMesh&, Istream&)",
                    is
                )   << "Base gradient scheme "
                    << baseGradSchemeName_
                    << " is not supported!"
                    << exit(FatalIOError);
            }

            tinterpScheme_ =
            tmp<surfaceInterpolationScheme<Type> >
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            );
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

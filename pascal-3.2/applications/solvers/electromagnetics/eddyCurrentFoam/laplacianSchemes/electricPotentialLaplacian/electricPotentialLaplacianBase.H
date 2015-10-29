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
    Foam::fv::electricPotentialLaplacianBase

Description
    TODO

SourceFiles
    electricPotentialLaplacianBase.C

\*---------------------------------------------------------------------------*/

#ifndef electricPotentialLaplacianBase_H
#define electricPotentialLaplacianBase_H

#include "word.H"
#include "Istream.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class electricPotentialLaplacianBase Declaration
\*---------------------------------------------------------------------------*/

class electricPotentialLaplacianBase
{

public:

    //- Enumeration defining the valid
    //  base laplacian schemes
    enum scheme
    {
        GAUSS,
        UNDEFINED
    };

private:

    // Private data

        //- Base laplacian scheme
        scheme scheme_;

        //- Names of valid base laplacian schemes
        static const NamedEnum<scheme, 1> schemeNames_;

public:

    // Constructors

        //- Construct from scheme
        electricPotentialLaplacianBase(const scheme scheme)
        {
            setBaseScheme(scheme);
        }

        //- Construct from scheme name
        electricPotentialLaplacianBase(const word schemeName)
        {
            setBaseScheme(schemeName);
        }

    // Destructor

        virtual ~electricPotentialLaplacianBase()
        {}

    // Member Functions

        //- Set base laplacian scheme
        void setBaseScheme(const scheme scheme)
        {
            scheme_ = scheme;
        }

        //- Set base laplacian scheme
        void setBaseScheme(const word schemeName)
        {
            if (schemeName == "")
            {
                FatalErrorIn
                (
                    "electricPotentialLaplacianBase::setBaseScheme(const word)"
                )   << "Base laplacian scheme not specified!" << endl
                    << "Valid base laplacian schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            typename NamedEnum<scheme, 1>::const_iterator
                cIter = schemeNames_.find(schemeName);

            if (cIter == schemeNames_.end())
            {
                FatalErrorIn
                (
                    "electricPotentialLaplacianBase::setBaseScheme(const word)"
                )   << "Base laplacian scheme " << schemeName
                    << " is not supported!" << endl
                    << "Valid base laplacian schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }
            else
            {
                setBaseScheme(schemeNames_[schemeName]);
            }
        }

        //- Set base laplacian scheme
        void setBaseScheme(Istream& is)
        {
            setBaseScheme(word(is));
        }

        //- Check base laplacian scheme
        void checkBaseScheme() const
        {
            if (scheme_ == UNDEFINED)
            {
                FatalErrorIn
                (
                    "electricPotentialLaplacianBase::checkBaseScheme()"
                )   << "Base laplacian scheme not specified!" << endl
                    << "Valid base laplacian schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }
        }

        //- Return base laplacian scheme
        scheme baseScheme() const
        {
            return scheme_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

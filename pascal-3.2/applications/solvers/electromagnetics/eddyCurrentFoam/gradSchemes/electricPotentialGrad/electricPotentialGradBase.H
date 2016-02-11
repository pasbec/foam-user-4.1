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
    Foam::fv::electricPotentialGradBase

Description
    TODO

SourceFiles
    electricPotentialGradBase.C

\*---------------------------------------------------------------------------*/

#ifndef electricPotentialGradBase_H
#define electricPotentialGradBase_H

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
               Class electricPotentialGradBase Declaration
\*---------------------------------------------------------------------------*/

class electricPotentialGradBase
{

public:

    //- Enumeration defining the valid
    //  base gradient schemes
    enum scheme
    {
        GAUSS,
        LEASTSQUARES,
        UNDEFINED
    };

private:

    // Private data

        //- Base gradient scheme
        scheme scheme_;

        //- Names of valid base gradient schemes
        static const NamedEnum<scheme, 2> schemeNames_;

public:

    // Constructors

        //- Construct from scheme
        electricPotentialGradBase(const scheme scheme)
        {
            setBaseScheme(scheme);
        }

        //- Construct from scheme name
        electricPotentialGradBase(const word schemeName)
        {
            setBaseScheme(schemeName);
        }

    // Destructor

        virtual ~electricPotentialGradBase()
        {}

    // Member Functions

        //- Set base gradient scheme
        void setBaseScheme(const scheme scheme)
        {
            scheme_ = scheme;
        }

        //- Set base gradient scheme
        void setBaseScheme(const word schemeName)
        {
            if (schemeName == "")
            {
                FatalErrorIn
                (
                    "electricPotentialGradBase::setBaseScheme(const word)"
                )   << "Base gradient scheme not specified!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            typename NamedEnum<scheme, 2>::const_iterator
                cIter = schemeNames_.find(schemeName);

            if (cIter == schemeNames_.end())
            {
                FatalErrorIn
                (
                    "electricPotentialGradBase::setBaseScheme(const word)"
                )   << "Base gradient scheme " << schemeName
                    << " is not supported!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }
            else
            {
                setBaseScheme(schemeNames_[schemeName]);
            }
        }

        //- Set base gradient scheme
        void setBaseScheme(Istream& is)
        {
            setBaseScheme(word(is));
        }

        //- Check base gradient scheme
        void checkBaseScheme()
        {
            if (scheme_ == UNDEFINED)
            {
                FatalErrorIn
                (
                    "electricPotentialGradBase::checkBaseScheme()"
                )   << "Base gradient scheme not specified!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }
        }

        //- Return base gradient scheme
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

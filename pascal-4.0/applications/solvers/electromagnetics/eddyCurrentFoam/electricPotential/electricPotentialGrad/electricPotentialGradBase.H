/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    ...

SourceFiles
    electricPotentialGradBase.C

Author
    Pascal Beckstein

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
    enum BaseScheme
    {
        GAUSS,
        LEASTSQUARES,
        UNDEFINED
    };


private:

    // Private data

        //- Base gradient scheme
        BaseScheme scheme_;

        //- Names of valid base gradient schemes
        static const NamedEnum<BaseScheme, 2> schemeNames_;


protected:

    // Protected Member Functions

        //- Set base gradient scheme
        void setBaseScheme(const word& schemeName)
        {
            if (schemeName == "")
            {
                FatalErrorIn
                (
                    "electricPotentialGradBase::setBaseScheme(const word&)"
                )   << "Base gradient scheme not specified!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }

            typename NamedEnum<BaseScheme, 2>::const_iterator
                cIter = schemeNames_.find(schemeName);

            if (cIter == schemeNames_.end())
            {
                FatalErrorIn
                (
                    "electricPotentialGradBase::setBaseScheme(const word&)"
                )   << "Base gradient scheme " << schemeName
                    << " is not supported!" << endl
                    << "Valid base gradient schemes are :" << endl
                    << schemeNames_.sortedToc()
                    << exit(FatalIOError);
            }
            else
            {
                scheme_ = schemeNames_[schemeName];
            }
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
        BaseScheme baseScheme() const
        {
            return scheme_;
        }

public:

    // Constructors

        //- Construct from scheme
        electricPotentialGradBase(BaseScheme scheme)
        :
            scheme_(scheme)
        {}

        //- Construct from scheme name
        electricPotentialGradBase(const word& schemeName)
        {
            setBaseScheme(schemeName);
        }

        //- Construct from stream
        electricPotentialGradBase(Istream& is)
        {
            setBaseScheme(word(is));
        }

    // Destructor

        virtual ~electricPotentialGradBase()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

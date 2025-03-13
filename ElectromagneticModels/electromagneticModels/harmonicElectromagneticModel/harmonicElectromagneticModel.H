/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::harmonicElectromagneticModel

Description
    Derived-class from electromagneticModel.

See also
    Foam::electromagneticModel

SourceFiles
    harmonicElectromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef harmonicElectromagneticModel_H
#define harmonicElectromagneticModel_H

#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class harmonicElectromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class harmonicElectromagneticModel
:
    virtual public electromagneticModel
{

    //- Electric potential [V]
        //- real
        volScalarField& PotEre_;
        //- imaginary
        volScalarField& PotEim_;

    //- Current density correction [A/m^2]
        //real
        volVectorField& deltaJre_;
        //imaginary
        volVectorField& deltaJim_;

    //- Velocity contribution [T*m/s]
        //real
        volVectorField& deltaUxBre_;
        //imaginary
        volVectorField& deltaUxBim_;

    //- Current density [A/m^2]
        //- real
        volVectorField& Jre_;
        // Not the actual current, but for keeping reference value
        volVectorField& JreReference_;

        //- imaginary
        volVectorField& Jim_;
        // Not the actual current, but for keeping reference value
        volVectorField& JimReference_;

    //- Magnetic field [T]
        //- real
        volVectorField& Bre_;

        //- imaginary
        volVectorField& Bim_;

    //- isComplex_ switch for solvers
        // true for harmonic models
        // fields have imaginary components
        const bool isComplex_ = true;

public:


    //- Runtime type information
    TypeName("harmonicElectromagneticModel");

    // Constructors

        //- Construct from components
        harmonicElectromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        harmonicElectromagneticModel(const harmonicElectromagneticModel&) = delete;

    //- Destructor
    virtual ~harmonicElectromagneticModel();


    // Member Functions

        //- Electric potential [V]
        virtual volScalarField& PotE(bool imaginary=false) override;

        //- Current density [A/m^2]
        virtual volVectorField& J(bool imaginary=false) override;
        virtual volVectorField& Jref(bool imaginary=false) override;

        //- Magnetic field [T]
        virtual volVectorField& B(bool imaginary=false) override;

        //- Update velocity contribution
        virtual void updateDeltaU(volVectorField& Udiff) override;

        //- Current density correction [A/m^2]
        virtual volVectorField& deltaJ(bool imaginary=false) override;

        //- For determining if model fields has imaginary parts
        virtual bool isComplex() const override;

        //const-access

        //- Electric potential [V]
        const virtual volScalarField& PotE(bool imaginary=false) const override;

        //- Current density [A/m^2]
        const virtual volVectorField& J(bool imaginary=false) const override;
        const virtual volVectorField& Jref(bool imaginary=false) const override;

        //- Magnetic field [T]
        const virtual volVectorField& B(bool imaginary=false) const override;

        //- Velocity contribution [T*m/s]
        const virtual volVectorField& deltaUxB(bool imaginary=false) const override;

        //- Current density correction [A/m^2]
        const virtual volVectorField& deltaJ(bool imaginary=false) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Complied within electromagneticModels.C
#ifdef NoRepository
    #include "harmonicElectromagneticModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::transientElectromagneticModel

Description
    Base-class for fluid thermodynamic properties.

See also
    Foam::electromagneticModel

SourceFiles
    transientElectromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef transientElectromagneticModel_H
#define transientElectromagneticModel_H

#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class transientElectromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class transientElectromagneticModel
:
    virtual public electromagneticModel
{

    //- Electric potential [V]
    volScalarField& PotE_;

    //- Current density correction [A/m^2]
    volVectorField& deltaJre_;

    //- Current density [A/m^2]
    //const char* JName_ = "Jre";
    volVectorField& Jre_;

    //- Magnetic field [T]
    //const char* BName_ = "Bre";
    volVectorField& Bre_;

    //- isComplex_ switch for solvers
    // false for transient models
    const bool isComplex_ = false;

public:


    //- Runtime type information
    TypeName("transientElectromagneticModel");


    // Constructors

        //- Construct from components
        transientElectromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );


    // Selectors
/*
        //- Standard selection based on fvMesh
        static autoPtr<transientElectromagneticModel> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );
*/
        //- Disallow default bitwise copy construction
        transientElectromagneticModel(const transientElectromagneticModel&) = delete;


    //- Destructor
    virtual ~transientElectromagneticModel();


    // Member Functions

        //- Electric potential [V]
        virtual volScalarField& PotE(bool imaginary=false) override;

        //- Current density [A/m^2]
        virtual volVectorField& J(bool imaginary=false) override;

        //- Magnetic field [T]
        virtual volVectorField& B(bool imaginary=false) override;

        //- Current density correction [A/m^2]
        virtual volVectorField& deltaJ(bool imaginary=false) override;

        //- For determining if model fields has imaginary parts
        virtual bool isComplex() const override;

        //const-access

        //- Electric potential [V]
        const virtual volScalarField& PotE(bool imaginary=false) const override;

        //- Current density [A/m^2]
        const virtual volVectorField& J(bool imaginary=false) const override;

        //- Magnetic field [T]
        const virtual volVectorField& B(bool imaginary=false) const override;

        //- Current density correction [A/m^2]
        const virtual volVectorField& deltaJ(bool imaginary=false) const override;

        virtual word getCoupledPotentialName(const word Jname) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Complied within electromagneticModels.C
#ifdef NoRepository
    #include "transientElectromagneticModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Base-class for fluid thermodynamic properties.

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

    //- Current density [A/m^2]
        //- real
        //const char* JreName_ = "Jre";
        volVectorField& Jre_;

        //- imaginary
        //const char* JimName_ = "Jim";
        volVectorField& Jim_;

    //- Magnetic field [T]
        //- real
        //const char* BreName_ = "Bre";
        volVectorField& Bre_;

        //- imaginary
        //const char* BimName_ = "Bim";
        volVectorField& Bim_;

    //- isComplex_ switch for solvers
        // true for harmonic models
        // fields have imaginary components
        const bool isComplex_ = true;

public:


    //- Runtime type information
    TypeName("harmonicElectromagneticModel");


    //- Declare run-time constructor selection table
/*
    declareRunTimeNewSelectionTable
    (
        autoPtr,
        harmonicElectromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );
*/

    // Constructors

        //- Construct from components
        harmonicElectromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );


    // Selectors
/*
        //- Standard selection based on fvMesh
        static autoPtr<harmonicElectromagneticModel> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );
*/
        //- Disallow default bitwise copy construction
        harmonicElectromagneticModel(const harmonicElectromagneticModel&) = delete;

    //- Destructor
    virtual ~harmonicElectromagneticModel();


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

        virtual const word getCoupledPotentialName(const word Jname) const override;
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

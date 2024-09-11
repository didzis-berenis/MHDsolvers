/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::electromagneticModel

Description
    Base-class for transientElectromagneticModel and harmonicElectromagneticModel.

    The electromagneticModel class is the pure virtual interface, plus static data and
    functions, including the selection table. It is default-constructable and
    holds no non-static data.

    The derived levels should be complete enough to act as the base of a
    electromagnetic instantiation.

SourceFiles
    electromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticModel_H
#define electromagneticModel_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class electromagneticModel Declaration
\*---------------------------------------------------------------------------*/

//This is an abstract class
class electromagneticModel
:
    public IOdictionary
{
    // Private Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        // Switch for updating electromagnetics

        bool correctElectromagnetics_ = false;

protected:

    // Protected data

        //- Phase-name
        const word& phaseName_;

        //- Electrical conductivity [S/m]
        const dimensionedScalar sigmaConst_;

        // Fields

        // Lorentz force
        volVectorField& JxB_;

        // Joule heating
        volScalarField& JJsigma_;

        //- Electrical conductivity [S/m]
        volScalarField& sigma_;

        //- Inverse Electric conductivity [m/S]
        volScalarField& sigmaInv_;

        //- Velocity difference
        //volVectorField& deltaU_;


    // Protected Static Member Functions

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volScalarField& lookupOrConstructScalar
        (
            const fvMesh& mesh,
            const char* name,
            writeOption wo = IOobject::AUTO_WRITE
        );

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volVectorField& lookupOrConstructVector
        (
            const fvMesh& mesh,
            const char* name,
            writeOption wo = IOobject::AUTO_WRITE
        );

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volScalarField& lookupOrConstructScalar
        (
            const fvMesh& mesh,
            const char* name,
            dimensionedScalar value,
            readOption ro = IOobject::NO_READ,
            writeOption wo = IOobject::AUTO_WRITE//Change to AUTO_WRITE for debugging
        );

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volVectorField& lookupOrConstructVector
        (
            const fvMesh& mesh,
            const char* name,
            dimensionedVector value,
            readOption ro = IOobject::NO_READ,
            writeOption wo = IOobject::AUTO_WRITE//Change to AUTO_WRITE for debugging
        );

        //- Lookup the named field, or construct it as a copy of the given field
        static volVectorField& lookupOrConstructVector
        (
            const fvMesh& mesh,
            const char* name,
            const volVectorField& field,
            readOption ro = IOobject::NO_READ,
            writeOption wo = IOobject::AUTO_WRITE//Change to AUTO_WRITE for debugging
        );

        //- Construct named field as MUST-READ if it is not found
        static void constructScalar
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Construct named field as MUST-READ if it is not found
        static void constructVector
        (
            const fvMesh& mesh,
            const char* name
        );

        static typeIOobject<IOdictionary> readModelDict
        (
            const objectRegistry& obr,
            const word& group,
            bool registerObject = false
        );


public:

    //- Runtime type information
    TypeName("electromagneticModel");

    //- Declare run-time constructor selection table
    // define electromagneticModel::fvMeshConstructorTable
    declareRunTimeSelectionTable
    (
        autoPtr,
        electromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    //- Generic New for each of the related electromagnetics model
    static inline autoPtr<electromagneticModel> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );

    // Static Member Functions

        //- Lookup the thermo associated with the given patch field
        static const electromagneticModel& lookupElectromagnetic(const fvPatchScalarField& pf);

    // Constructors

        //- Construct from components
        electromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName=word::null
        );

        //- Disallow default bitwise copy construction
        electromagneticModel(const electromagneticModel&) = delete;


    //- Destructor
    virtual ~electromagneticModel()
    {}

    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Set correctElectromagnetics switch to true
        virtual void setCorrectElectromagnetics();

        //- Set correctElectromagnetics switch to false
        virtual void setCorrected();

        //- Get correctElectromagnetics switch value
        virtual bool correctElectromagnetics() const;

        //- calculate and return J correction difference
        virtual void solve();

        //- calculate and return J correction difference
        virtual void findPotE(bool imaginary=false);

        //- calculate and return J correction difference
        virtual void findDeltaJ(bool imaginary=false);

        //- Inverse Electric conductivity [m/S]
        virtual volScalarField& sigmaInv();

        //- Electric conductivity [S/m]
        virtual volScalarField& sigma();

        //- assign closest internal field values to sigma boundaries
        //virtual void patchSigmaBoundaries(volScalarField& sigma);

        //- return field from object registry
        virtual volScalarField& getScalarFieldRef(const char* name);

        //- return field from object registry
        virtual volVectorField& getVectorFieldRef(const char* name);

        //- Predict the electromagnetic model coefficients if possible
        //  without solving electromagnetic model equations
        virtual void predict();

        //- Solve the electromagnetic model equations
        //- Update properties
        virtual void correct();

        virtual void updateDeltaU(volVectorField& Udiff) = 0;

        //Template functions for derived classes

        //- For determining if model fields has imaginary parts
        virtual bool isComplex() const = 0;

        //- Electric potential [V]
        virtual volScalarField& PotE(bool imaginary=false) = 0;

        //- Current density [A/m^2]
        virtual volVectorField& J(bool imaginary=false) = 0;

        //- Magnetic field [T]
        virtual volVectorField& B(bool imaginary=false) = 0;

        //- Current density correction [A/m^2]
        virtual volVectorField& deltaJ(bool imaginary=false) = 0;

        virtual word getCoupledCurrentName(const word Jname) const = 0;

        //Read-only access

        //- Properties dictionary
        //virtual const IOdictionary& properties() const;

        //- Return const access to the mesh
        //virtual const fvMesh& mesh() const;

        //- return field from object registry
        const virtual volScalarField& getScalarField(const char* name) const;

        //- return field from object registry
        const virtual volVectorField& getVectorField(const char* name) const;

        //- Electric conductivity [S/m]
        virtual const volScalarField& sigma() const;

        //- Electric conductivity for patch [S/m]
        virtual tmp<scalarField> sigma(const label patchi) const;

        //- Inverse Electric conductivity for patch [S/m]
        //virtual tmp<scalarField> sigmaInv(const label patchi) const;

        //- Inverse Electric conductivity [m/S]
        virtual const volScalarField& sigmaInv() const;

        //- Electric conductivity [S/m]
        virtual const dimensionedScalar sigmaConst() const;

        //- Electric potential [V]
        virtual const volScalarField& PotE(bool imaginary=false) const = 0;

        //- Current density [A/m^2]
        virtual const volVectorField& J(bool imaginary=false) const = 0;

        //- Magnetic field [T]
        virtual const volVectorField& B(bool imaginary=false) const = 0;

        //- Magnetic field [T]
        virtual const volVectorField& deltaUxB(bool imaginary=false) const = 0;

        //- Current density correction [A/m^2]
        virtual const volVectorField& deltaJ(bool imaginary=false) const = 0;

        //- Lorentz force [N/m^3]
        const volVectorField& JxB;

        // Joule heating [W/m^3]
        const volScalarField& JJsigma;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "electromagneticModelNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

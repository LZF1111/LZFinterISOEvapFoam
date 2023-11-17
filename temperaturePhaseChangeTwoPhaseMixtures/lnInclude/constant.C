/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/
#include "Time.H"
#include "constant.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "twoPhaseMixtureEThermo.H"
#include "fvmSup.H"
#include "wallPolyPatch.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace temperaturePhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        temperaturePhaseChangeTwoPhaseMixture,
        constant,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::constant
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    temperaturePhaseChangeTwoPhaseMixture(mixture, mesh),
    coeffC_
    (
        "coeffC",
        dimless/dimTime/dimTemperature,
        optionalSubDict(type() + "Coeffs")
    ),
    coeffE_
    (
        "coeffE",
        dimless/dimTime/dimTemperature*dimLength,
        optionalSubDict(type() + "Coeffs")
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotAlphal() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
   const volScalarField& maggradalpha1 = mesh_.lookupObject<volScalarField>("maggradalpha1");
 // const volScalarField& alpha1 = mesh_.lookupObject<volScalarField>("alpha1");
    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);
 // volVectorField gradAlpha ( mixture_.alpha1()*fvc::grad(mixture_.alpha2()) - mixture_.alpha2()*fvc::grad(mixture_.alpha1()));
  //volVectorField gradAlpha1 (fvc::grad(mixture_.alpha1()));
  //volScalarField MagA ( mag(gradAlpha));
  //volScalarField MagA1 ( mag(gradAlpha1));
  /*volVectorField gradAlpha2
  (
        IOobject
        (
            "gradAlpha2 ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ, //若文件存在，则读取，否则
            IOobject::AUTO_WRITE       // 保存该变量
        ),
       
        gradAlpha
    );
   volVectorField gradAlpha3
  (
        IOobject
        (
            "gradAlpha3 ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ, //若文件存在，则读取，否则
            IOobject::AUTO_WRITE       // 保存该变量
        ),
       fvc::grad(mixture_.alpha1())
    );*/
 
  //volScalarField MagA ( mag(gradAlpha1));
//Info<< "magAlphagradient has caculated"<< MagA <<endl  ;

//Info<< "magAlphagradient has caculated"<<endl  ;

//Info<< "magAlphagradient has caculated"<< MagA <<endl  ;


    
    
    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*max(TSat - T, T0),
       -coeffE_*maggradalpha1*mixture_.rho1()*max(T - TSat, T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDot() const
{

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T"); 
  
    const volScalarField& maggradalpha1 = mesh_.lookupObject<volScalarField>("maggradalpha1");   
       
    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    const dimensionedScalar T0(dimTemperature, Zero);
    // const volScalarField& alpha1 = mesh_.lookupObject<volScalarField>("alpha1");
    volScalarField mDotE
    (
        "mDotE", coeffE_*maggradalpha1*mixture_.rho1()*limitedAlpha1*max(T - TSat, T0)
    );
    volScalarField mDotC
    (
        "mDotC", coeffC_*mixture_.rho2()*limitedAlpha2*max(TSat - T, T0)
    );

    if (mesh_.time().outputTime())
    {
        mDotC.write();
        mDotE.write();
    }

    return Pair<tmp<volScalarField>>
    (
        tmp<volScalarField>(new volScalarField(mDotC)),
        tmp<volScalarField>(new volScalarField(-mDotE))
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::mDotDeltaT() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );
   //  const volScalarField& alpha1 = mesh_.lookupObject<volScalarField>("alpha1");
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const volScalarField& maggradalpha1 = mesh_.lookupObject<volScalarField>("maggradalpha1");
    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    return Pair<tmp<volScalarField>>
    (
        coeffC_*mixture_.rho2()*limitedAlpha2*pos(TSat - T),
        coeffE_*maggradalpha1*mixture_.rho1()*limitedAlpha1*pos(T - TSat)
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::TSource() const
{
      //const volScalarField& alpha1 = mesh_.lookupObject<volScalarField>("alpha1");
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const volScalarField& maggradalpha1 = mesh_.lookupObject<volScalarField>("maggradalpha1");
    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();

    dimensionedScalar L = mixture_.Hf2() - mixture_.Hf1();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField Vcoeff
    (
        coeffE_*maggradalpha1*mixture_.rho1()*limitedAlpha1*L*pos(T - TSat)
    );
    const volScalarField Ccoeff
    (
        coeffC_*mixture_.rho2()*limitedAlpha2*L*pos(TSat - T)
    );

    TSource =
        fvm::Sp(Vcoeff, T) - Vcoeff*TSat
      + fvm::Sp(Ccoeff, T) - Ccoeff*TSat;

    return tTSource;
}



void Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::correct()
{
}


bool Foam::temperaturePhaseChangeTwoPhaseMixtures::constant::read()
{
    if (temperaturePhaseChangeTwoPhaseMixture::read())
    {
        subDict(type() + "Coeffs").readEntry("coeffC", coeffC_);
        subDict(type() + "Coeffs").readEntry("coeffE", coeffE_);

        return true;
    }

    return false;
}


// ************************************************************************* //

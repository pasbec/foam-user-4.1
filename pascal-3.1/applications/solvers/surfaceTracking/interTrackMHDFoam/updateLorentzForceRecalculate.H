
    // Write basic fields to have them in paraview while waiting for update
    U.write();
    p.write();

    // Write fields
    FL.write();
    pB.write();


//     // Basic security check for system call
//     if
//     (
//         (
//             lorentzForceUpdateUtility.find("rm ") > 0
//         )
//         ||
//         (
//             lorentzForceUpdateUtility.find("rmdir ") > 0
//         )
//         ||
//         (
//             lorentzForceUpdateUtility.find("-delete") > 0
//         )
//     )
//     {
//         FatalError
//             << "Update utility contains dangerous patterns (rm, rmdir, -delete)!"
//             << abort(FatalError);
//     }

    // Create system call string and inform
    string sysCallString = lorentzForceUpdateUtility + " " + runTime.timeName();
    Info << "Calling:" << " "
         << "\"" << sysCallString << "\"" << endl;

    // Update Lorentz force via system call
    int sysCallReturn = std::system(sysCallString.c_str());

    if(sysCallReturn != 0)
    {
        FatalError
            << "System call operation failed!"
            << abort(FatalError);
    }


    // Re-read Lorentz force FL
    Info << "Rereading field FL from time folder" << endl;

    volVectorField* FLuPtr = NULL;

    FLuPtr = new volVectorField
    (
        IOobject
        (
            "FL",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    volVectorField& FLu = *FLuPtr;

    FL == FLu; // == Copies boundary values, too

    delete FLuPtr; FLuPtr = NULL;


    // Re-read Magnetic pressure pB
    Info << "Rereading field pB from time folder" << endl;

    volScalarField* pBuPtr = NULL;

    pBuPtr = new volScalarField
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    volScalarField& pBu = *pBuPtr;

    pB == pBu; // == Copies boundary values, too

    delete pBuPtr; pBuPtr = NULL;



//     // ===========================================================================//
//     // Constants
//     // ===========================================================================//
//
//     dimensionedScalar mu0
//     (
//         "mu0",
//         dimensionSet(1, 1, -2, 0, 0, -2, 0),
//         4.0*mathematicalConstant::pi*1e-07
//     );
//
//     // ===========================================================================//
//     // FIXME [High]: This coil model is only valid for rgs2D!
//     // ===========================================================================//
//
//     scalar coilDx    = 0.005;
//     scalar coilDy    = 0.005;
//     scalar coilDz    = 0.1; // Grid z-size
//     scalar coilI     = 1000.0;
//     scalar coilF     = 1000.0;
//
//     List<vector> segmentCentres
//     (
//         IStringStream
//         (
//             "("
//                 "(-0.065 0.0096 0)"
//                 "(-0.065 0.0226 0)"
//                 "(-0.065 0.0356 0)"
//                 "( 0.065 0.0096 0)"
//                 "( 0.065 0.0226 0)"
//                 "( 0.065 0.0356 0)"
//             ")"
//         )()
//     );
//
//     List<vector> jS, jSvolMu0ByPi;
//     jS.setSize(segmentCentres.size());
//     jSvolMu0ByPi.setSize(segmentCentres.size());
//
//     forAll(segmentCentres, segmentI)
//     {
//         scalar coilIRMS = Foam::sqrt(2.0) * coilI;
//         scalar coilA = coilDx * coilDy;
//
//         vector jSdirSegmentI = vector(0,0,sign(segmentCentres[segmentI].component(0)));
//
//         jS[segmentI] = coilIRMS * jSdirSegmentI / coilA;
//
//         jSvolMu0ByPi[segmentI] = 1.0e-7 * coilIRMS * jSdirSegmentI * coilDz;
//     }
//
//     // ======================================================================//
//     // FIXME [High]: Material properties for liquid silicon
//     // ======================================================================//
//
//     dimensionedScalar sigma
//     (
//         "sigma",
//         dimensionSet(-1,-3,3,0,0,2,0),
//         scalar(1.2e+6)
//     );
//
//     // ===========================================================================//
//     // Frequency, Omega and alpha
//     // ===========================================================================//
//
//     dimensionedScalar freq
//     (
//         "freq",
//         dimensionSet(0,0,-1,0,0,0,0),
//         scalar(coilF)
//     );
//
//     dimensionedScalar omega = mathematicalConstant::twoPi*freq;
//
//     dimensionedScalar alpha = mu0 * sigma * omega;
//
//     // ===========================================================================//
//     // Biot-Savart law to get boundary conditions from source current density
//     // ===========================================================================//
//
//     volVectorField AReS
//     (
//         IOobject
//         (
//             "AReS",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh,
//         dimensionedVector
//         (
//             "zero",
//             dimMass/dimLength/dimCurrent/pow(dimTime,2),
//             vector::zero
//         ),
//         fixedValueFvPatchVectorField::typeName
//     );
//
//     forAll(mesh.boundary(), patchI)
//     {
//         const vectorField& faceCentres = mesh.boundaryMesh()[patchI].faceCentres();
//
//         if (!mesh.boundary()[patchI].coupled())
//         {
//             const fvPatch& patch = mesh.boundary()[patchI];
//
//             vector dr,AReSfaceI = vector(0,0,0);
//
//             forAll(patch, faceI)
//             {
//                 AReSfaceI = vector(0,0,0);
//
//                 forAll(segmentCentres, segmentI)
//                 {
//                     dr = faceCentres[faceI] - segmentCentres[segmentI];
//
//                     AReSfaceI += jSvolMu0ByPi[segmentI] / mag(dr);
//
// //                     AReSfaceI += ( jSvolMu0ByPi[segmentI] ^ dr) / pow(mag(dr),3);
//                 }
//
//                 AReS.boundaryField()[patchI][faceI] = AReSfaceI;
//             }
//         }
//     }
//
//     // ======================================================================//
//     // Get external boundary part for ARe and AIm
//     // ======================================================================//
//
//     forAll(mesh.boundary(), patchI)
//     {
//         if (!mesh.boundary()[patchI].coupled())
//         {
//             ARe.boundaryField()[patchI] == AReS.boundaryField()[patchI];
//             AIm.boundaryField()[patchI] == vectorField(AIm.boundaryField()[patchI].size(), vector::zero);
//         }
//     }
//
//     // ===========================================================================//
//     // Iterative solution to calculate full fields ARe and AIm
//     // ===========================================================================//
//
//     int nBcorr = 11;
//     for (int Bcorr=0; Bcorr<nBcorr; Bcorr++)
//     {
//         // ==================================================================//
//         // Block matrix preparations
//         // ==================================================================//
//
//         // Prepare block system
//         BlockLduMatrix<vector6> AEqn(mesh);
//
//         // Reference block diagonal and set it to zero
//         Field<tensor6>& AEqnD = AEqn.diag().asSquare(); AEqnD = tensor6::zero;
//
//         // Reference linear off-diagonal and set it to zero
//         Field<vector6>& AEqnU = AEqn.upper().asLinear(); AEqnU = vector6::zero;
//         Field<vector6>& AEqnL = AEqn.lower().asLinear(); AEqnL = vector6::zero;
//
//         // Set up source
//         vector6Field AEqnS(mesh.nCells(), vector6::zero);
//
//         forAll(mesh.boundary(), patchI)
//         {
//             if (!mesh.boundary()[patchI].coupled())
//             {
//                 vector6Field& patch = A.boundaryField()[patchI];
//
//                 patch.component(0) = ARe.boundaryField()[patchI].component(0);
//                 patch.component(1) = ARe.boundaryField()[patchI].component(1);
//                 patch.component(2) = ARe.boundaryField()[patchI].component(2);
//                 patch.component(3) = AIm.boundaryField()[patchI].component(0);
//                 patch.component(4) = AIm.boundaryField()[patchI].component(1);
//                 patch.component(5) = AIm.boundaryField()[patchI].component(2);
//             }
//         }
//
//         // Transfer the coupled interface list for processor/cyclic/etc. boundaries
//         AEqn.interfaces() = A.boundaryField().blockInterfaces();
//
//         // ==================================================================//
//         // Real part
//         // ==================================================================//
//
//         // Discretize real part
//         tmp<fvVectorMatrix> AReEqn
//         (
//             - fvm::laplacian(ARe)
//             - alpha * AIm
//         );
//
//         // ==================================================================//
//         // Imaginary part
//         // ==================================================================//
//
//         tmp<fvVectorMatrix> AImEqn
//         (
//             - fvm::laplacian(AIm)
//             + alpha * ARe
//         );
//
//         // ==================================================================//
//         // Assemble segregated equations in block structure
//         // ==================================================================//
//
//         // Inset equations into block Matrix
//         blockMatrixTools::insertEquation(0, AReEqn(), AEqn, A.internalField(), AEqnS);
//         blockMatrixTools::insertEquation(3, AImEqn(), AEqn, A.internalField(), AEqnS);
//
//         // ==================================================================//
//         // Assemble point-implicit coupled equation part in block structure
//         // ==================================================================//
//
//         // Add off-diagonal terms and remove them from block source
//         forAll(AEqnD, i)
//         {
//
//             AEqnD[i](0,3) = -alpha.value() * mesh.V()[i];
//             AEqnD[i](3,0) = +alpha.value() * mesh.V()[i];
//             AEqnD[i](1,4) = -alpha.value() * mesh.V()[i];
//             AEqnD[i](4,1) = +alpha.value() * mesh.V()[i];
//             AEqnD[i](2,5) = -alpha.value() * mesh.V()[i];
//             AEqnD[i](5,2) = +alpha.value() * mesh.V()[i];
//
//             AEqnS[i][0] -= alpha.value() * A.internalField()[i][3] * mesh.V()[i];
//             AEqnS[i][3] += alpha.value() * A.internalField()[i][0] * mesh.V()[i];
//             AEqnS[i][1] -= alpha.value() * A.internalField()[i][4] * mesh.V()[i];
//             AEqnS[i][4] += alpha.value() * A.internalField()[i][1] * mesh.V()[i];
//             AEqnS[i][2] -= alpha.value() * A.internalField()[i][5] * mesh.V()[i];
//             AEqnS[i][5] += alpha.value() * A.internalField()[i][2] * mesh.V()[i];
//         }
//
//         // ==================================================================//
//         // Solve the block matrix
//         // ==================================================================//
//
//         BlockSolverPerformance<vector6> solverPerf =
//             BlockLduSolver<vector6>::New
//             (
//                 A.name(),
//                 AEqn,
//                 mesh.solutionDict().solver(A.name())
//             )->solve(A,AEqnS);
//
//         solverPerf.print();
//
//         // ==================================================================//
//         // Retrieve the solution and update for next iteration
//         // ==================================================================//
//
//         blockMatrixTools::retrieveSolution(0, ARe.internalField(), A);
//         blockMatrixTools::retrieveSolution(3, AIm.internalField(), A);
//
//         ARe.correctBoundaryConditions();
//         ARe.correctBoundaryConditions();
//
//         AReEqn.clear();
//         AImEqn.clear();
//
// //         // =======================================================================//
// //         // Segregated approach for inner solution
// //         // =======================================================================//
// //
// //         Vector<label>& geoP = const_cast<Vector<label>&>(mesh.geometricD());
// //         Vector<label>& solP = const_cast<Vector<label>&>(mesh.solutionD());
// //
// //         Vector<label> geoPback = geoP;
// //         Vector<label> solPback = solP;
// //
// //         for (int i=0; i<3;i++)
// //         {
// //             geoP.component(i) = 1;
// //             solP.component(i) = 1;
// //         }
// //
// //         int nBsegCorr = 10;
// //         for (int BsegCorr=0; BsegCorr<nBsegCorr; BsegCorr++)
// //         {
// //             fvVectorMatrix AReEqn
// //             (
// //                 - fvm::laplacian(ARe)
// //                 - mu0 * sigma * omega * AIm
// //             );
// //
// //             fvVectorMatrix AImEqn
// //             (
// //                 - fvm::laplacian(AIm)
// //                 + mu0 * sigma * omega * ARe
// //             );
// //
// //             AReEqn.solve();
// //             AImEqn.solve();
// //
// //             ARe.correctBoundaryConditions();
// //             AIm.correctBoundaryConditions();
// //         }
// //
// //         geoP = geoPback;
// //         solP = solPback;
//
//         // ==================================================================//
//         // Calculate current density
//         // ==================================================================//
//
//         jRe ==  sigma * omega * AIm;
//         jIm == -sigma * omega * ARe;
//
//         BRe == fvc::curl(ARe);
//         BIm == fvc::curl(AIm);
//
// //         jRe == (1.0/mu0) * fvc::curl(BRe);
// //         jIm == (1.0/mu0) * fvc::curl(BIm);
//
//         // ==================================================================//
//         // Correct boundary values based on Biot-Savart law
//         // ==================================================================//
//
//         vectorField jReVolMu0ByPi = 1.0e-7 * jRe.internalField() * mesh.V();
//         vectorField jImVolMu0ByPi = 1.0e-7 * jIm.internalField() * mesh.V();
//
//         const vectorField& cellCentres = mesh.C();
//
//         forAll(mesh.boundary(), patchI)
//         {
//             if (!mesh.boundary()[patchI].coupled())
//             {
//                 const fvPatch& patch = mesh.boundary()[patchI];
//
//                 const vectorField& faceCentres = mesh.boundaryMesh()[patchI].faceCentres();
//
//                 vector dr,AReFaceI, AImFaceI;
// //                 scalar drMagPow3;
//
//                 ARe.boundaryField()[patchI] == vectorField(ARe.boundaryField()[patchI].size(), vector::zero);
// //                 AIm.boundaryField()[patchI] == vectorField(AIm.boundaryField()[patchI].size(), vector::zero);
//
//                 forAll(patch, faceI)
//                 {
//                     AReFaceI = vector(0,0,0);
//                     AImFaceI = vector(0,0,0);
//
//                     forAll(cellCentres, cellI)
//                     {
//                         dr = faceCentres[faceI] - cellCentres[cellI];
//
// //                         drMagPow3 = 1.0 / pow(mag(dr),3);
//
// //                         AReFaceI += (jReVolMu0ByPi[cellI] ^ dr) * drMagPow3;
// //                         AImFaceI += (jImVolMu0ByPi[cellI] ^ dr) * drMagPow3;
//
//                         AReFaceI += jReVolMu0ByPi[cellI] / mag(dr);
// //                         AImFaceI += jImVolMu0ByPi[cellI] / mag(dr);
//                     }
//
//                     ARe.boundaryField()[patchI][faceI] = AReFaceI;
// //                     AIm.boundaryField()[patchI][faceI] = AImFaceI;
//                 }
//
//                 ARe.boundaryField()[patchI] == AReS.boundaryField()[patchI]
//                                              + ARe.boundaryField()[patchI];
//             }
//         }
//     }
//
//     // ===========================================================================//
//     // Calculate FL and pB
//     // ===========================================================================//
//
// //     FL == (1.0/rho) * (0.5 * ((jRe ^ ARe) + (jIm ^ AIm)));
// //     pB == (0.5/mu0) * (0.5 * ((ARe & ARe) + (AIm & AIm)));

//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010-2012 Alessandro Tasora
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be 
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//      
//   Demo code about   
//   
//     - collisions and contacts 
//
//       (This is just a possible method of integration
//       of Chrono::Engine + Irrlicht: many others 
//       are possible.)
//     
//	 CHRONO   
//   ------
//   Multibody dinamics engine 
//   
// ------------------------------------------------ 
//             www.deltaknowledge.com
// ------------------------------------------------ 
///////////////////////////////////////////////////
    
  
 
#include "physics/ChApidll.h" 
#include "physics/ChSystem.h"
#include "irrlicht_interface/ChBodySceneNode.h"
#include "irrlicht_interface/ChBodySceneNodeTools.h" 
#include "irrlicht_interface/ChIrrAppInterface.h"
#include "core/ChRealtimeStep.h"
#include "lcp/ChLcpIterativeMINRES.h" // test
#include "physics/ChMaterialSurface.h"
#include <cstring>

 
// Use the namespace of Chrono

using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;


const double rho = 1000;

void Calc_Hydrodynamics_Force(ChVector<float> & force3, ChVector<float> & forceLoc,
		ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {

	// ***** calculation of force
	ChVector<> freeSurfaceNormal = -mphysicalSystem.Get_G_acc();
	float g = freeSurfaceNormal.Length();
	freeSurfaceNormal.Normalize();

	ChVector<> bodyCtr = mrigidBody->GetPos();
	ChVector<> dist3 = freeSurfaceLocation - bodyCtr;
	double dist = dist3.Dot(freeSurfaceNormal);

	float rad = mrigidBody->GetCollisionModel()->GetSafeMargin(); //this only works for sphere

	force3 = ChVector<>(0,0,0);
	forceLoc = bodyCtr;
	if (dist > rad) { //outside water
		return;
	} else if (dist < -rad) {
		double V = 4.0 / 3 * CH_C_PI * pow(rad, 3);
		force3 = V * rho * g;
		forceLoc = bodyCtr;
	} else {
		double h = rad - dist;
		double a = sqrt(rad * rad - dist * dist);
		double V = CH_C_PI * h * h / 3 * (3 * rad - h);
		force3 = V * rho * g;
		double distFromCenter = 3.0 / 4 * pow(2 * rad - h, 2) / (3 * rad - h); 	// http://mathworld.wolfram.com/SphericalCap.html -->
																				// Harris and Stocker 1998, p. 107 (Harris, J. W. and Stocker,
																				// H. "Spherical Segment (Spherical Cap)." §4.8.4 in Handbook of
																				// Mathematics and Computational Science. New York: Springer-Verlag, p. 107, 1998.)
		forceLoc = distFromCenter * (-freeSurfaceNormal);
	}
	force3 = ChVector<>(0,0,0);
}

void add_hydronynamic_force(ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {
	// ***** insertion of force
	ChSharedPtr<ChForce> hydroForce;
//	std::string forceTag("hydrodynamics");
	char forceTag[] = "hydrodynamics";
	hydroForce = mrigidBody->SearchForce(forceTag);
	if (hydroForce.IsNull()) {
		printf("hydro create\n");

		hydroForce = ChSharedPtr<ChForce>(new ChForce);
		mrigidBody->AddForce(hydroForce);
		// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
		hydroForce->SetName(forceTag);
	}
	ChVector<float> force3;
	ChVector<float> forceLoc;

	Calc_Hydrodynamics_Force(force3, forceLoc, mrigidBody, mphysicalSystem, freeSurfaceLocation);

	hydroForce->SetVrelpoint(forceLoc);
	hydroForce->SetMforce(force3.Length());
	force3.Normalize();
	hydroForce->SetDir(force3);
}

void update_hydronynamic_force(ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {
	// ***** insertion of force
	ChSharedPtr<ChForce> hydroForce;
//	std::string forceTag("hydrodynamics");
	char forceTag[] = "hydrodynamics";
	hydroForce = mrigidBody->SearchForce(forceTag);
	if (!hydroForce.IsNull()) {
		printf("hydro update\n");
		ChVector<float> force3;
		ChVector<float> forceLoc;
		Calc_Hydrodynamics_Force(force3, forceLoc, mrigidBody, mphysicalSystem, freeSurfaceLocation);
		force3 = ChVector<float>(0,0,0);

		hydroForce->SetVrelpoint(forceLoc);
		hydroForce->SetMforce(force3.Length());
		force3.Normalize();
		hydroForce->SetDir(force3);
	}
}

void create_some_falling_items(ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
{

	ChBodySceneNode* mrigidBody; 

	ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
	mmaterial->SetFriction(0.4f);
	mmaterial->SetCompliance (0.00005f);
	mmaterial->SetComplianceT(0.00005f);
	mmaterial->SetDampingF(0.2);
	mmaterial->SetCohesion(0.05);

	// Create a bunch of ChronoENGINE rigid bodies (spheres and
	// boxes) which will fall..
	// Bodies are Irrlicht nodes of the special class ChBodySceneNode, 
	// which encapsulates ChBody items).  
	
	video::ITexture* cubeMap   = driver->getTexture("../data/cubetexture_borders.png");
	video::ITexture* sphereMap = driver->getTexture("../data/bluwhite.png");
		
// Arman: Blocks deactivated
//	for (int ai = 0; ai < 1; ai++)  // N. of walls
//	{
//		for (int bi = 0; bi < 10; bi++)  // N. of vert. bricks
//		{
//			for (int ui = 0; ui < 15; ui++)  // N. of hor. bricks
//			{
//				mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
//													&mphysicalSystem, msceneManager,
//													0.8,
//													ChVector<>(-8+ui*4.0+2*(bi%2),  1.0+bi*2.0, ai*9),
//													ChQuaternion<>(1,0,0,0),
//													ChVector<>(3.96,2,4) );
//				mrigidBody->GetBody()->SetMaterialSurface(mmaterial);
//				mrigidBody->setMaterialTexture(0,	cubeMap);
//				mrigidBody->addShadowVolumeSceneNode();
//			}
//		}
//	}


	// Jenga tower
	/*
	for (int bi = 0; bi < 12; bi+=2) 
	{ 
		ChBodySceneNode* mrigidBody1;
		mrigidBody1 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											0.2,
											ChVector<>(-5, 1.0+bi*2.0,  0),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(2,2, 14) );
		mrigidBody1->GetBody()->SetMaterialSurface(mmaterial);
		mrigidBody1->setMaterialTexture(0,	cubeMap);

		ChBodySceneNode* mrigidBody2;
		mrigidBody2 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											0.2,
											ChVector<>( 5, 1.0+bi*2.0,  0),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(2,2, 14) );
		mrigidBody2->GetBody()->SetMaterialSurface(mmaterial);
		mrigidBody2->setMaterialTexture(0,	cubeMap);

		ChBodySceneNode* mrigidBody3;
		mrigidBody3 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											0.2,
											ChVector<>(0, 3.0+bi*2.0,  5),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(14,2, 2) );
		mrigidBody3->GetBody()->SetMaterialSurface(mmaterial);
		mrigidBody3->setMaterialTexture(0,	cubeMap);

		ChBodySceneNode* mrigidBody4;
		mrigidBody4 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											0.2,
											ChVector<>(0, 3.0+bi*2.0,  -5),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(14,2, 2) );
		mrigidBody4->GetBody()->SetMaterialSurface(mmaterial);
		mrigidBody4->setMaterialTexture(0,	cubeMap);
	}
	*/

	// Create the floor using
	// fixed rigid body of 'box' type:

	mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											200.0,
											ChVector<>(0,-2,0),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(550,4,550) );
	mrigidBody->GetBody()->SetBodyFixed(true);
	mrigidBody->GetBody()->SetMaterialSurface(mmaterial);


	// Create a ball that will collide with wall
	double mradius = 4;
	double density = 1.01;
	double mmass = (4./3.)*CH_C_PI*pow(mradius,3)*density; 
	GetLog() << "Ball mass = " << mmass << "\n";
	double minert = (2./5.)* mmass * pow(mradius,2);

	mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
										&mphysicalSystem, msceneManager,
										mmass, // mass
										ChVector<>(0, 10, -8), // pos
										mradius, // radius
										20,  // hslices, for rendering
										15); // vslices, for rendering

	// set moment of inertia (more realistic than default 1,1,1).
	mrigidBody->GetBody()->SetInertiaXX(ChVector<>(minert,minert,minert));
	mrigidBody->GetBody()->SetPos_dt(ChVector<>(0,0,0));
	mrigidBody->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	mrigidBody->GetBody()->GetMaterialSurface()->SetCompliance(0.0);
	mrigidBody->GetBody()->GetMaterialSurface()->SetComplianceT(0.0);
	mrigidBody->GetBody()->GetMaterialSurface()->SetDampingF(0.2);

	// Some aesthetics for 3d view..
	mrigidBody->addShadowVolumeSceneNode();
	mrigidBody->setMaterialTexture(0,	sphereMap);

//	void add_hydronynamic_force(ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {


//	add_hydronynamic_force(mrigidBody->GetBody().get_ptr(), mphysicalSystem, ChVector<>(0, 5, -8));

}


 

 
 
int main(int argc, char* argv[])
{ 

	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed.
	DLL_CreateGlobals();

	// Create a ChronoENGINE physical system
	ChSystem mphysicalSystem; 

	// Create the Irrlicht visualization (open the Irrlicht device, 
	// bind a simple user interface, etc. etc.)
	ChIrrAppInterface application(&mphysicalSystem, L"Bricks test",core::dimension2d<u32>(800,600),false, true); 

 
	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo  (application.GetDevice());
	ChIrrWizard::add_typical_Sky   (application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 590,  400);
	ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15,14,-30), core::vector3df(0,5,0)); 

	// 
	// HERE YOU CREATE THE MECHANICAL SYSTEM OF CHRONO... 
	// 

 
	// Create all the rigid bodies.
	create_some_falling_items(mphysicalSystem, application.GetSceneManager(), application.GetVideoDriver());
  
  
	// Prepare the physical system for the simulation 

	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR_MULTITHREAD);

	mphysicalSystem.SetUseSleeping(false);

	mphysicalSystem.SetMaxPenetrationRecoverySpeed(1.6); // used by Anitescu stepper only
	mphysicalSystem.SetIterLCPmaxItersSpeed(40);
	mphysicalSystem.SetIterLCPmaxItersStab(20); // unuseful for Anitescu, only Tasora uses this
	mphysicalSystem.SetIterLCPwarmStarting(true);
	mphysicalSystem.SetParallelThreadNumber(4);

	//
	// THE SOFT-REAL-TIME CYCLE
	//
 
	application.SetStepManage(true);
	application.SetTimestep(0.02);

	while(application.GetDevice()->run())
	{
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));

		ChIrrTools::drawGrid(application.GetVideoDriver(), 5,5, 20,20, 
			ChCoordsys<>(ChVector<>(0,0.2,0),Q_from_AngAxis(CH_C_PI/2,VECT_X)), video::SColor(50,90,90,150),true);

		application.DrawAll();

		application.DoStep();
 
		application.GetVideoDriver()->endScene();  

		//************

//		for(int i=0; i<mphysicalSystem.Get_bodylist()->size(); i++){
//			update_hydronynamic_force(mphysicalSystem.Get_bodylist()->at(i), mphysicalSystem, ChVector<>(0, 5, -8));
//
//		}

//		std::vector<ChBody*>::iterator ibody = mphysicalSystem.Get_bodylist()->begin();
//		while (ibody != mphysicalSystem.Get_bodylist()->end()) {
//			update_hydronynamic_force(*ibody, mphysicalSystem, ChVector<>(0, 5, -8));
//			ibody++;
////			printf("body pos %f %f %f\n", (*ibody)->coord.pos.x, (*ibody)->coord.pos.y, (*ibody)->coord.pos.z);
//		}

	}
	return 0;
}
  

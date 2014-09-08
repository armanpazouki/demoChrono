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
    
  
 
//#include "physics/ChApidll.h"
#include "physics/ChSystem.h"
#include "unit_IRRLICHT/ChBodySceneNode.h"
#include "unit_IRRLICHT/ChBodySceneNodeTools.h"
#include "unit_IRRLICHT/ChIrrAppInterface.h"
#include "core/ChRealtimeStep.h"
#include "lcp/ChLcpIterativeMINRES.h" // test
#include "physics/ChMaterialSurface.h"
#include "physics/ChController.h"
#include "collision/ChCModelBulletBody.h"
#include "physics/ChContactContainer.h"
#include "core/ChTimer.h"
#include <cstring>
#include <fstream>
#include <map>

 
// Use the namespace of Chrono

using namespace chrono;
using namespace collision;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;
using namespace std;



const double rhoF = 1000;
const double rhoR = 917;
const double mu_Viscosity = 1;//.1;
const ChVector<> surfaceLoc = ChVector<>(0, .09, -.08);

//******************* ship and sphere stuff
double mradius = .02;
int numLayers = 5;

ChBodySceneNode* shipPtr;
const double shipVelocity = .0027;//.27;//1; //arman modify
double shipInitialPosZ = 0;
const double timePause = .1; //arman modify
double ship_width = .20;
double box_X = ship_width, box_Y = .50, box_Z = .02;
double collisionEvelop = .06 * mradius;
//**********************************

void Calc_Hydrodynamics_Forces(ChVector<> & F_Hydro, ChVector<> & forceLoc, ChVector<> & T_Drag,
		ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation) {
	F_Hydro = ChVector<>(0,0,0);
	forceLoc = ChVector<>(0,0,0);
	T_Drag = ChVector<>(0,0,0);


	// ***** calculation of force
	ChVector<> freeSurfaceNormal = -mphysicalSystem.Get_G_acc();
	double g = freeSurfaceNormal.Length();
	freeSurfaceNormal.Normalize();

	ChVector<> bodyCtr = mrigidBody->GetPos();
	ChVector<> dist3 = bodyCtr - freeSurfaceLocation;
	double dist = dist3.Dot(freeSurfaceNormal); // distance of the sphere center from the fluid surface

	double rad = mrigidBody->GetCollisionModel()->GetSafeMargin(); //this only works for sphere

	//****************** Buoyancy Force
	ChVector<> F_Buoyancy = ChVector<>(0,0,0);
	forceLoc = bodyCtr;
	if (dist < -rad) {
		double V = 4.0 / 3 * CH_C_PI * pow(rad, 3);
		F_Buoyancy = V * rhoF * g * freeSurfaceNormal;
		forceLoc = bodyCtr;
	} else if (dist < rad) {
		double h = rad - dist;
		double V = CH_C_PI * h * h / 3 * (3 * rad - h);
		F_Buoyancy = V * rhoF * g * freeSurfaceNormal;
		double distFromCenter = 3.0 / 4 * pow(2 * rad - h, 2) / (3 * rad - h); 	// http://mathworld.wolfram.com/SphericalCap.html -->
																				// Harris and Stocker 1998, p. 107 (Harris, J. W. and Stocker,
																				// H. "Spherical Segment (Spherical Cap)." §4.8.4 in Handbook of
																				// Mathematics and Computational Science. New York: Springer-Verlag, p. 107, 1998.)
		forceLoc = bodyCtr + distFromCenter * (-freeSurfaceNormal);
	}
	// "dist > rad" --> outside of water
	//****************** Drag Force and Torque
	double Cd = 0.4;
	ChVector<> vel = mrigidBody->GetPos_dt();
	ChVector<> F_Drag = ChVector<>(0,0,0);
	if (dist < rad) {
		double A_ref = 0.5 * CH_C_PI * rad * (rad - dist);
		double multDrag = 1;
		if (mphysicalSystem.GetChTime() < timePause) {
			multDrag = 1;
		} else {
			multDrag = 1;
		}
		F_Drag = multDrag * (-6.0 * CH_C_PI * mu_Viscosity * rad * vel
					-0.5 * rhoF * Cd * vel.Length() * vel);
		T_Drag = -8.0 * CH_C_PI * mu_Viscosity * pow(rad, 3) * mrigidBody->GetWvel_par(); // in parent, i.e. absoute, reference frame.
	}
	//****************** Total Force
	F_Hydro = F_Buoyancy + F_Drag; // it is assumed that F_Drag is applied at the buoyancy center
}

void create_hydronynamic_force(ChBody* mrigidBody, ChSystem& mphysicalSystem, const chrono::ChVector<>& freeSurfaceLocation, bool createForce) {
	// ***** insertion of force
	ChSharedPtr<ChForce> hydroForce;
	ChSharedPtr<ChForce> hydroTorque;

//	std::string forceTag("hydrodynamics_force");
	char forceTag[] = "hydrodynamics_force";
	char torqueTag[] = "hydrodynamics_torque";
	hydroForce = mrigidBody->SearchForce(forceTag);
	hydroTorque = mrigidBody->SearchForce(torqueTag);

	//********** create force if needed **********
	if (createForce) {
		if (hydroForce.IsNull()) {
			hydroForce = ChSharedPtr<ChForce>(new ChForce);
			hydroForce->SetMode(FTYPE_FORCE); // no need for this. It is the default option.
			mrigidBody->AddForce(hydroForce);
			// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
			hydroForce->SetName(forceTag);
		}
		if (hydroTorque.IsNull()) {
			hydroTorque = ChSharedPtr<ChForce>(new ChForce);
			hydroTorque->SetMode(FTYPE_TORQUE);
			mrigidBody->AddForce(hydroTorque);
			// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
			hydroTorque->SetName(torqueTag);
		}
	}
	//********** update force magnitude **********
	if (!hydroForce.IsNull() || !hydroTorque.IsNull()) {
		ChVector<> F_Hydro;
		ChVector<> forceLoc;
		ChVector<> T_Drag;

		Calc_Hydrodynamics_Forces(F_Hydro, forceLoc, T_Drag, mrigidBody, mphysicalSystem, freeSurfaceLocation);

		hydroForce->SetVpoint(forceLoc);
		hydroForce->SetMforce(F_Hydro.Length());
		F_Hydro.Normalize();
		hydroForce->SetDir(F_Hydro);

		hydroTorque->SetMforce(T_Drag.Length());
		T_Drag.Normalize();
		hydroTorque->SetDir(T_Drag);
	}
}

void calc_ship_contact_forces(ChSystem& mphysicalSystem, ChVector<> & mForce, ChVector<> & mTorque) {
	mForce = ChVector<>(0,0,0);
	mTorque = ChVector<>(0,0,0);
	ChContactContainer* container  = (ChContactContainer *) mphysicalSystem.GetContactContainer();
//	std::map<ChBody*, ChVector<> > m_forces;
//	std::map<ChBody*, ChVector<> > m_torques;
//	ChVector<> mForce;
//	ChVector<> mTorque;

	std::list<ChContact*> m_list = container->GetContactList();
	for (std::list<ChContact *>::iterator it=m_list.begin(); it != m_list.end(); ++it){
	  ChVector<> force_contactFrame = (*it)->GetContactForce();
	  ChVector<> force_abs = *((*it)->GetContactPlane()) * force_contactFrame;
	  ChModelBulletBody * model_A = (ChModelBulletBody *) (*it)->GetModelA();
	  ChModelBulletBody * model_B = (ChModelBulletBody *) (*it)->GetModelB();
//		  if ((model_A->GetBody() != shipPtr->GetBody()) && (model_B->GetBody() != shipPtr->GetBody())) {
//			  continue;
//		  }
	  ChBody * body_A = model_A->GetBody();
	  ChBody * body_B = model_B->GetBody();

	  if (body_A == (ChBody*)shipPtr->GetBody().get_ptr()) {
		  mForce -= force_abs;

		  ChVector<> point_on_A = (*it)->GetContactP1();
		  mTorque -= (point_on_A - body_A->GetPos()) % force_abs;
//		  ChVector<> local_point_on_A = ChTransform<>::TransformParentToLocal(point_on_A, body_A->GetPos(), body_A->GetRot());
//		  mTorque += local_point_on_A % force;
	  } else if (body_B == (ChBody*)shipPtr->GetBody().get_ptr()) {
		  mForce += force_abs;

		  ChVector<> point_on_B = (*it)->GetContactP2();
		  mTorque += (point_on_B - body_B->GetPos()) % force_abs;
//		  ChVector<> local_point_on_B = ChTransform<>::TransformParentToLocal(point_on_B, body_B->GetPos(), body_B->GetRot());
//		  mTorque -= local_point_on_B % force;
	  }
	}
}


void create_ice_particles(ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
{

	ChBodySceneNode* mrigidBody;
	ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
	mmaterial->SetFriction(0.4f);

		//mmaterial->SetCompliance (0.00005f);
		//mmaterial->SetComplianceT(0.00005f);
		//mmaterial->SetDampingF(0.2);
		//mmaterial->SetCohesion(0.05);

	// Create a bunch of ChronoENGINE rigid bodies (spheres and
	// boxes) which will fall..
	// Bodies are Irrlicht nodes of the special class ChBodySceneNode, 
	// which encapsulates ChBody items).  
	
	video::ITexture* cubeMap   = driver->getTexture("../data/cubetexture_borders.png");
	video::ITexture* sphereMap = driver->getTexture("../data/bluwhite.png");

	// Create the floor using
	// fixed rigid body of 'box' type:

	//*** create bed
	ChBodySceneNode* earthPtr;
	earthPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>(0,-.40,0),
											ChQuaternion<>(1,0,0,0), 
											ChVector<>(5.50,.04,5.50) );
	earthPtr->GetBody()->SetBodyFixed(true);
	earthPtr->GetBody()->SetMaterialSurface(mmaterial);
	earthPtr->GetBody()->SetCollide(false);

	//*****
	ChVector<> boxMin = ChVector<>(-.04, .05, -.12);
	ChVector<> boxMax = ChVector<>(-.04 + .80, .05, -.12 + 1.10);
	//**************** add walls
	//*** side wall 1
	double wall_width = .60;
	double ship_height = .005;
	ChBodySceneNode* wallPtr;
	wallPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>(boxMin.x - ship_height/2, 0, .50),
											ChQuaternion<>(1,0,0,0),
											ChVector<>(ship_height,wall_width,1.5) );
	wallPtr->GetBody()->SetBodyFixed(true);
	wallPtr->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	wallPtr->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

	//*** side wall 2
	wallPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>(boxMax.x + ship_height/2, 0, .50),
											ChQuaternion<>(1,0,0,0),
											ChVector<>(ship_height,wall_width,1.5) );
	wallPtr->GetBody()->SetBodyFixed(true);
	wallPtr->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	wallPtr->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

	//*** end wall
	wallPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>((boxMin.x + boxMax.x)/2, 0, boxMax.z),
											ChQuaternion<>(1,0,0,0),
											ChVector<>(boxMax.x - boxMin.x, wall_width, ship_height) );
	wallPtr->GetBody()->SetBodyFixed(true);
	wallPtr->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	wallPtr->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

	//*** beginning walls
	double hole_width = 1.2 * ship_width;
	double small_wall_Length = 0.5 * (boxMax.x - boxMin.x - hole_width);
	wallPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>(boxMin.x + .5 * small_wall_Length, 0, boxMin.z - ship_height/2),
											ChQuaternion<>(1,0,0,0),
											ChVector<>( small_wall_Length, wall_width, ship_height) );
	wallPtr->GetBody()->SetBodyFixed(true);
	wallPtr->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	wallPtr->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

	wallPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
											&mphysicalSystem, msceneManager,
											rhoR,
											ChVector<>(boxMax.x - .5 * small_wall_Length, 0, boxMin.z - ship_height/2),
											ChQuaternion<>(1,0,0,0),
											ChVector<>( small_wall_Length, wall_width, ship_height) );
	wallPtr->GetBody()->SetBodyFixed(true);
	wallPtr->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
	wallPtr->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default


	//**************** sphere prob
	// note: dimensions are in cm
	double spacing = 2 * mradius*1.1;
	int numColX = (boxMax.x - boxMin.x - spacing) / spacing;
	int numColZ = (boxMax.z - boxMin.z - spacing) / spacing;

	for (int j = 0; j < numLayers; j++) {
		for (int i = 0; i < numColX; i++) {
			for (int k = 0; k < numColZ; k++) {
				// Create a ball that will collide with wall
				double mmass = (4./3.)*CH_C_PI*pow(mradius,3)*rhoR;
				double minert = (2./5.)* mmass * pow(mradius,2);

				mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
													&mphysicalSystem, msceneManager,
													mmass, // mass
													ChVector<>(i * spacing, j * spacing, k * spacing) + (boxMin + .5 * ChVector<>(spacing,spacing,spacing)) , // pos
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
				mrigidBody->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

				// Some aesthetics for 3d view..
				//mrigidBody->addShadowVolumeSceneNode();
				mrigidBody->setMaterialTexture(0,	sphereMap);

				create_hydronynamic_force(mrigidBody->GetBody().get_ptr(), mphysicalSystem, surfaceLoc, true);
			}
		}

	}

//	// Create a ball that will collide with wall
//	double mmass = (4./3.)*CH_C_PI*pow(mradius,3)*rhoR;
//	double minert = (2./5.)* mmass * pow(mradius,2);
//
//	mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
//										&mphysicalSystem, msceneManager,
//										mmass, // mass
//										ChVector<>(.5 * (boxMax.x + boxMin.x),  .09, shipInitialPosZ + 2*mradius),
//										mradius, // radius
//										20,  // hslices, for rendering
//										15); // vslices, for rendering
//
//	// set moment of inertia (more realistic than default 1,1,1).
//	mrigidBody->GetBody()->SetInertiaXX(ChVector<>(minert,minert,minert));
//	mrigidBody->GetBody()->SetPos_dt(ChVector<>(0,0,0));
//	mrigidBody->GetBody()->GetMaterialSurface()->SetFriction(0.4f);
//	mrigidBody->GetBody()->GetMaterialSurface()->SetCompliance(0.0);
//	mrigidBody->GetBody()->GetMaterialSurface()->SetComplianceT(0.0);
//	mrigidBody->GetBody()->GetMaterialSurface()->SetDampingF(0.2);
//
//	// Some aesthetics for 3d view..
//	mrigidBody->addShadowVolumeSceneNode();
//	mrigidBody->setMaterialTexture(0,	sphereMap);
//
//	create_hydronynamic_force(mrigidBody->GetBody().get_ptr(), mphysicalSystem, surfaceLoc, true);

	//*** create ship
	double boxMass = rhoR * box_X * box_Y * box_Z;
	double bI1 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Y, 2));
	double bI2 = 1.0 / 12 * boxMass * (pow(box_Y, 2) + pow(box_Z, 2));
	double bI3 = 1.0 / 12 * boxMass * (pow(box_X, 2) + pow(box_Z, 2));

	shipInitialPosZ = boxMin.z - box_Z;
	shipPtr = (ChBodySceneNode*)addChBodySceneNode_easyBox(
										&mphysicalSystem, msceneManager,
										boxMass,
										ChVector<>(.5 * (boxMax.x + boxMin.x),  .09, shipInitialPosZ),
										ChQuaternion<>(1,0,0,0),
										ChVector<>(box_X, box_Y, box_Z) );
	shipPtr->GetBody()->SetMass(boxMass);
	shipPtr->GetBody()->SetInertiaXX(ChVector<>(bI2, bI3, bI1));
	shipPtr->GetBody()->SetMaterialSurface(mmaterial);
	shipPtr->setMaterialTexture(0,	cubeMap);
	wallPtr->GetBody()->GetCollisionModel()->SetEnvelope(collisionEvelop); //envelop is .03 by default

	shipPtr->addShadowVolumeSceneNode();
	shipPtr->GetBody()->SetPos_dt(ChVector<>(0,0,0));
	char forceTag[] = "pulling_force";
	ChSharedPtr<ChForce> pullingForce = ChSharedPtr<ChForce>(new ChForce);
	pullingForce->SetMode(FTYPE_FORCE); // no need for this. It is the default option.
	shipPtr->GetBody()->AddForce(pullingForce);
	// ** or: hydroForce = ChSharedPtr<ChForce>(new ChForce());
	pullingForce->SetName(forceTag);
	pullingForce->SetVpoint(shipPtr->GetBody()->GetPos());
	pullingForce->SetMforce(0);
	pullingForce->SetDir(ChVector<>(1,0,0));

	//***** prismatic constraint between ship and bed
//	ChSharedPtr<ChLinkLockPlanePlane> shipConstraint(new ChLinkLockPlanePlane);
//	shipConstraint->Initialize(shipPtr->GetBody(), earthPtr->GetBody(),
//			ChCoordsys<>(ChVector<>(30,  9, -25) , Q_from_AngAxis(CH_C_PI/2, VECT_X))
//			);
	ChSharedPtr<ChLinkLockPrismatic> shipConstraint(new ChLinkLockPrismatic);
	shipConstraint->Initialize(shipPtr->GetBody(), earthPtr->GetBody(),
			ChCoordsys<>(ChVector<>(.30,  .09, -.25) , QUNIT)
			);
	mphysicalSystem.AddLink(shipConstraint);
}

void MoveShip(ChSystem& mphysicalSystem) {
	static bool onCall = false;
//	if (!onCall) {
//		onCall = true;
		shipPtr->GetBody()->SetPos_dt(ChVector<>(0,0,shipVelocity));
//	}
//    ChSharedPtr<ChControllerPID> my_controllerPID(new ChControllerPID);
//    my_controllerPID->P = 1.0e9;
//    my_controllerPID->D = 1.0e8;
//    my_controllerPID->I = 1.0e8;
//
//    double forcePID_X = my_controllerPID->Get_Out(shipPtr->GetBody()->GetPos().z - shipInitialPosZ - shipVelocity * (mphysicalSystem.GetChTime() - timePause), mphysicalSystem.GetChTime());
//    char forceTag[] = "pulling_force";
//	ChSharedPtr<ChForce> pullingForce = shipPtr->GetBody()->SearchForce(forceTag);
//	pullingForce->SetMforce(forcePID_X);
//	pullingForce->SetDir(ChVector<>(0,0,-1));
}

void FixShip(ChSystem& mphysicalSystem) {
	shipPtr->GetBody()->SetBodyFixed(true);
}
 
int main(int argc, char* argv[])
{ 
	ChTimer<double> myTimer;
	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed.
//	DLL_CreateGlobals();

	// Create a ChronoENGINE physical system
	ChSystem mphysicalSystem; 

	// Create the Irrlicht visualization (open the Irrlicht device,
	// bind a simple user interface, etc. etc.)
	ChIrrAppInterface application(&mphysicalSystem, L"Bricks test",core::dimension2d<u32>(800,600),false, true);
	fstream outForceData("forceData.txt", ios::out);
// 1*********
  // *** Irrlicht stuff, deactivated
	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	ChIrrWizard::add_typical_Logo  (application.GetDevice());
	ChIrrWizard::add_typical_Sky   (application.GetDevice());
	ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(.7f, 2.2f, -.9f), core::vector3df(-3.0f, 8.0f, 6.0f), 59,  40);
	ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-.15,.4,-.40), core::vector3df(0,.05,0));
// 2*********

	// 
	// HERE YOU CREATE THE MECHANICAL SYSTEM OF CHRONO... 
	// 

 
	// Create all the rigid bodies.
	create_ice_particles(mphysicalSystem, application.GetSceneManager(), application.GetVideoDriver());
  
  
	// Prepare the physical system for the simulation 

	mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR);
	mphysicalSystem.SetUseSleeping(false);
	mphysicalSystem.SetMaxPenetrationRecoverySpeed(.003); // used by Anitescu stepper only
	mphysicalSystem.SetIterLCPmaxItersSpeed(500);
	//mphysicalSystem.SetIterLCPmaxItersStab(20); // unuseful for Anitescu, only Tasora uses this
	//mphysicalSystem.SetIterLCPwarmStarting(true);
	//mphysicalSystem.SetParallelThreadNumber(2);

	mphysicalSystem.SetTol(0);
	mphysicalSystem.SetTolSpeeds(.01 * shipVelocity);
	//
	// THE SOFT-REAL-TIME CYCLE
	//
 
	application.SetStepManage(true);
	application.SetTimestep(.001);  //Arman modify
//std::cout<<"reay to simulate"<<std::endl;

	outForceData << "time, forceX, forceY, forceZ, forceMag, pressureX, pressureY, pressureZ, pressureMag, shipVelocity, energy, timePerStep.## numSpheres" << mphysicalSystem.Get_bodylist()->end() - mphysicalSystem.Get_bodylist()->begin()
			<< " pauseTime: " << timePause<< " setVelocity: "<< shipVelocity << endl;
	while(application.GetDevice()->run() && mphysicalSystem.GetChTime() < 20) //arman modify
	{
		myTimer.start();
// 1********* irrlicht initialization 8888
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));
		ChIrrTools::drawGrid(application.GetVideoDriver(), .05,.05, 40,40,
			ChCoordsys<>(ChVector<>(0,-.30,0),Q_from_AngAxis(CH_C_PI/2,VECT_X)), video::SColor(50,90,90,150),true);
		application.DrawAll();
// 2*********
		application.DoStep();
//		std::cout<<"after step"<<std::endl;

// 1*********
		application.GetVideoDriver()->endScene();
// 2*********
		if (mphysicalSystem.GetChTime() > timePause) {
			MoveShip(mphysicalSystem);
			//shipPtr->GetBody()->SetPos_dt(ChVector<>(0,0,shipVelocity));
		} else {
			FixShip(mphysicalSystem);
		}
		//******************** ship force*********************
//		ChVector<> shipForce = shipPtr->GetBody()->Get_Xforce();
//		printf("force %f\n",shipForce.z);

		ChVector<> mForce;
		ChVector<> mTorque;
		calc_ship_contact_forces(mphysicalSystem, mForce, mTorque);
		ChVector<> icePressure = mForce / (numLayers * 2 * mradius * cos(CH_C_PI / 6)) / ship_width;


		myTimer.stop();
		//****************************************************
//		for(int i=0; i<mphysicalSystem.Get_bodylist()->size(); i++){
//			create_hydronynamic_force(mphysicalSystem.Get_bodylist()->at(i), mphysicalSystem, surfaceLoc, false);
//
//		}
		std::vector<ChBody*>::iterator ibody = mphysicalSystem.Get_bodylist()->begin();
		double energy = 0;
		while (ibody != mphysicalSystem.Get_bodylist()->end()) {
			create_hydronynamic_force(*ibody, mphysicalSystem, surfaceLoc, false);
			energy += pow((*ibody)->GetPos_dt().Length() , 2);
			ibody++;
		}
//		printf("time %f, force %f %f %f, shipVelocity %f, simulation time %f, energy %f\n", mphysicalSystem.GetChTime(), mForce.x, mForce.y, mForce.z, shipPtr->GetBody()->GetPos_dt().z, myTimer(), energy);
		outForceData << mphysicalSystem.GetChTime() << ", " << mForce.x << ", " << mForce.y << ", " << mForce.z << ", " <<
				mForce.Length() << ", " <<
				icePressure.x << ", " << icePressure.y << ", " << icePressure.z << ", " << icePressure.Length() << ", " <<
				shipPtr->GetBody()->GetPos_dt().z << ", " << energy << ", " << myTimer() << endl;

		printf("Time %f, energy %f\n", mphysicalSystem.GetChTime(), energy);
	}
	outForceData.close();
	return 0;
}
  

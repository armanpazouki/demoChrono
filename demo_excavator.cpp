#include "physics/ChApidll.h"
#include "physics/ChSystem.h"
#include "irrlicht_interface/ChBodySceneNode.h"
#include "irrlicht_interface/ChBodySceneNodeTools.h"
#include "irrlicht_interface/ChIrrAppInterface.h"
#include "core/ChRealtimeStep.h"
#include "physics/ChController.h"
#include <irrlicht.h>
#include <fstream>
#include "algorithm"
#include <ctime>

// Use the namespace of Chrono

using namespace chrono;
using namespace collision;
using namespace std;

// Use the main namespaces of Irrlicht
using namespace irr;
using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;

#define MIN(a,b) ((a)>(b)?(b):(a))
#define massScale 1
#define massScaleParticles 1

// returns 1.0f for positive floats, -1.0f for negative floats, 0.0f for zero
inline float fast_sign(float f) {
    if (((int&)f & 0x7FFFFFFF)==0) return 0.f; // test exponent & mantissa bits: is input zero?
    else {
        float r = 1.0f;
        (int&)r |= ((int&)f & 0x80000000); // mask sign bit in f, set it in r if necessary
        return r;
    }
}

// Static values valid through the entire program (bad
// programming practice, but enough for quick tests)

double STATIC_motor1 = 0;
double STATIC_motor2= 0;
double STATIC_motor3 = 0;
double STATIC_motor4 = 0;

//**************************************************************************************************************************
void PrintToFileForces(ChSharedPtr<ChLinkEngine> motor, double torque, const char* fileName, int step, double dT, int startWrite, int endWrite) {

        fstream outToFile;
        double time = dT * step;
        if (step == 1) {
                outToFile.open(fileName, ios::out);
                outToFile<<"step"<<"  "<<"Time"<<"  "<<"PID_torque"<<"  "<<
                        "reactForce_X"<<"  "<<"reactForce_Y"<<"  "<<"reactForce_Z"<<"  "<<
                        "global_reactForce_X"<<"  "<<"global_reactForce_X"<<"  "<<"global_reactForce_X"<<"  "<<
                        "react_torque_X"<<"  "<<"react_torque_Y"<<"  "<<"react_torque_Z"<<"  "<<
                        "global_react_torque_X"<<"  "<<"global_react_torque_X"<<"  "<<"global_react_torque_X"<<"  "<<endl;
        } else {
                outToFile.open(fileName, ios::app);
        }

        if (step > 2) {
                outToFile<<step<<"  "<<time<<"  "<<torque<<"  "<<
                        motor->Get_react_force().x<<"  "<<motor->Get_react_force().y<<"  "<<motor->Get_react_force().z<<"  "<<
                        motor->GetForce_X()->Get_iforce()<<"  "<<motor->GetForce_Y()->Get_iforce()<<"  "<<motor->GetForce_Z()->Get_iforce()<<"  "<<
                        motor->Get_react_torque().x<<"  "<<motor->Get_react_torque().y<<"  "<<motor->Get_react_torque().z<<"  "<<
                        motor->GetForce_Rx()->Get_iforce()<<"  "<<motor->GetForce_Ry()->Get_iforce()<<"  "<<motor->GetForce_Rz()->Get_iforce()<<"  "<<endl;
                //GetC_torque().z is equal to the torque at the previous step. I am guessing (as it is claimed in the header file CHlinkMarkers.h) that
                //GetC_force() and GetC_torque() return the applied external force/torque on the constraint.
        }

        outToFile.close();
}
//**************************************************************************************************************************
void PrintToFile_BodyForces(ChSharedBodyPtr& mbody, const char* fileName, int step) {

        fstream outToFile;
        double time = mbody->GetSystem()->GetStep() * step;
        if (step == 1) {
                outToFile.open(fileName, ios::out);
                outToFile<<"step"<<"  "<<"Time"<<"  "<<
                        "Get_accumulated_force_X"<<"  "<<"Get_accumulated_force_Y"<<"  "<<"Get_accumulated_force_Z"<<"  "<<
                        "Get_accumulated_torque_X"<<"  "<<"Get_accumulated_torque_Y"<<"  "<<"Get_accumulated_torque_Z"<<"  "<<
                        "Get_Xforce_X"<<"  "<<"Get_Xforce_Y"<<"  "<<"Get_Xforce_Z"<<"  "<<
                        "Get_Xtorque_X"<<"  "<<"Get_Xtorque_Y"<<"  "<<"Get_Xtorque_Z"<<"  "<<
                        "Get_Scr_force_X"<<"  "<<"Get_Scr_force_Y"<<"  "<<"Get_Scr_force_Z"<<"  "<<
                        "Get_Scr_torque_X"<<"  "<<"Get_Scr_torque_Y"<<"  "<<"Get_Scr_torque_Z"<<"  "<<
                        endl;
        } else {
                outToFile.open(fileName, ios::app);
        }

        if (step > 2) {
                outToFile<<step<<"  "<<time<<"  "<<
                        mbody->Get_accumulated_force().x<<"  "<<mbody->Get_accumulated_force().y<<"  "<<mbody->Get_accumulated_force().z<<"  "<<
                        mbody->Get_accumulated_torque().x<<"  "<<mbody->Get_accumulated_torque().y<<"  "<<mbody->Get_accumulated_torque().z<<"  "<<
                        mbody->Get_Xforce().x<<"  "<<mbody->Get_Xforce().y<<"  "<<mbody->Get_Xforce().z<<"  "<<
                        mbody->Get_Xtorque().x<<"  "<<mbody->Get_Xtorque().y<<"  "<<mbody->Get_Xtorque().z<<"  "<<
                        mbody->Get_Scr_force().x<<"  "<<mbody->Get_Scr_force().y<<"  "<<mbody->Get_Scr_force().z<<"  "<<
                        mbody->Get_Scr_torque().x<<"  "<<mbody->Get_Scr_torque().y<<"  "<<mbody->Get_Scr_torque().z<<"  "<<
                        endl;
        }

        outToFile.close();
}
//**************************************************************************************************************************
// motion controll function which can be used to perform a prescribed motion without the use of inverse kinematics
int MotionControl(double* rotation_control, int stage, int num_stages, int step, double* mem_control, ChSharedPtr<ChLinkEngine>  my_motor1, ChSharedPtr<ChLinkEngine>  my_motor2, ChSharedPtr<ChLinkEngine>  my_motor3, ChSharedPtr<ChLinkEngine>  my_motor4) {
        int myswitch=0;
        int num_step=0;
        double tol=0;
        if (stage>=num_stages)
                {
                // number of steps and tolerance for the return motion
                num_step=rotation_control[5];
                tol=rotation_control[4];
        }
        else
        {
                num_step=(int)rotation_control[stage*6+5];
                tol=rotation_control[stage*6+4];
        }
                                if (step<=1)
                                {
                                        //get previous contoller position
                                        mem_control[0]=STATIC_motor1;
                                        mem_control[1]=STATIC_motor2;
                                        mem_control[2]=STATIC_motor3;
                                        mem_control[3]=STATIC_motor4;
                                }

                                if (step<=num_step )
                                {
                                        // if all stages have been processed -> move back to zero
                                        if (stage>=num_stages)
                                        {
                                        // update the controller
                                        STATIC_motor1=mem_control[0]*(cos(CH_C_PI/num_step*step)+1.0)/2.0;
                                        STATIC_motor2=mem_control[1]*(cos(CH_C_PI/num_step*step)+1.0)/2.0;
                                        STATIC_motor3=mem_control[2]*(cos(CH_C_PI/num_step*step)+1.0)/2.0;
                                        STATIC_motor4=mem_control[3]*(cos(CH_C_PI/num_step*step)+1.0)/2.0;
                                        }
                                        else
                                        {
                                        // update the controller
                                        STATIC_motor1=mem_control[0]+rotation_control[stage*6+0]*(-cos(CH_C_PI/num_step*step)+1)/2;
                                        STATIC_motor2=mem_control[1]+rotation_control[stage*6+1]*(-cos(CH_C_PI/num_step*step)+1)/2;
                                        STATIC_motor3=mem_control[2]+rotation_control[stage*6+2]*(-cos(CH_C_PI/num_step*step)+1)/2;
                                        STATIC_motor4=mem_control[3]+rotation_control[stage*6+3]*(-cos(CH_C_PI/num_step*step)+1)/2;
                                        }
                                }
                                // check for controller success
                                if (fabs(my_motor1->Get_mot_rot()-STATIC_motor1)<tol && fabs(my_motor2->Get_mot_rot()-STATIC_motor2)<tol && fabs(my_motor3->Get_mot_rot()-STATIC_motor3)<tol && fabs(my_motor4->Get_mot_rot()-STATIC_motor4)<tol && step>=num_step)
                                {
                                        // switch to next motion stage
                                        myswitch=1;
                                }
        return myswitch;
}

////**************************************************************************************************************************
//// create the sand
//void create_some_falling_items(ChIrrAppInterface* app, ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
//{
//      ChBodySceneNode* mrigidBody;
//
//
//      for (int bi = 0; bi < 2000; bi ++) {
//              int width=63.0 / 2.0;
//              int height=11.0 / 2.0;
//              double Rs = 6.0 / 3.5 * 2;       // scoop 'bounding sphere' rad. for approximate inertia
//              double Ms = 4.0 / 3.0 *CH_C_PI*Rs*Rs*Rs*0.0000026 * massScaleParticles;  // sphere mass (granite)
//              double Js = (2.0/5.0)*Ms*Rs*Rs; // scoop XX YY ZZ inertia, as sphere, just to test things
//
//              mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
//                                                                                      &mphysicalSystem, msceneManager,
//                                                                                      Ms,
//                                                                                      ChVector<>(((bi-bi%height)%(height*width))/height*(-2.0f)*Rs-400,20+(bi%height)*2.0f*Rs,-16+Rs+2.0f*Rs*(bi-bi%(height*width))/(height*width)),
//                                                                                      Rs);
//
//              mrigidBody->GetBody()->SetFriction(0.40);
//              //mrigidBody->GetBody()->SetRollingFriction( 0.001 );
//              mrigidBody->GetBody()->SetInertiaXX(ChVector<>(Js,Js,Js));
//              //mrigidBody->addShadowVolumeSceneNode();
//
//
//              video::ITexture* sphereMap = driver->getTexture("../data/rock.jpg");
//              mrigidBody->setMaterialTexture(0,sphereMap);
//      }
//}
//**************************************************************************************************************************
// create the sand
void create_some_falling_items(ChIrrAppInterface* app, ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
{
        ChBodySceneNode* mrigidBody;
        double scaling = 2;
        int numParticles = 2000 * scaling;//2000;
        double Rs = 6.0 / 3.5 * 2 / pow(scaling, 1.0/3);         // scoop 'bounding sphere' rad. for approximate inertia
        double Ms = 4.0 / 3.0 *CH_C_PI*Rs*Rs*Rs*0.0000026 * massScaleParticles;  // sphere mass (granite)
        double Js = (2.0/5.0)*Ms*Rs*Rs; // scoop XX YY ZZ inertia, as sphere, just to test things

        int width= 107 / Rs;//63.0 / 2.0;
        int depth = 45 / Rs;//13;
        int height = numParticles / (width * depth) + 1;

        //int height=11.0 / 2.0;
        //int depth = numParticles / (width * height) + 1;
                //printf("************************* i am adding particle ********************* %d \n", depth);

        int bi = 0;
        for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) {
                        for (int i = 0; i < width; i++) {
                                if (bi < numParticles) {
                                        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
                                                                                                                                &mphysicalSystem, msceneManager,
                                                                                                                                Ms,
                                                                                                                                ChVector<>(-400 - i * 2.0 * Rs, 20 + j * 2.0f * Rs, -16 + Rs + 2.0f * Rs * k),
                                                                                                                                Rs);

                                        mrigidBody->GetBody()->SetFriction(0.40);
                                        //mrigidBody->GetBody()->SetRollingFriction( 0.001 );
                                                                mrigidBody->GetBody()->SetMass(Ms);
                                        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(Js,Js,Js));
                                        //mrigidBody->addShadowVolumeSceneNode();


                                        video::ITexture* sphereMap = driver->getTexture("../data/rock.jpg");
                                        mrigidBody->setMaterialTexture(0,sphereMap);
                                }
                                bi ++;
                        }
                }
        }
}
//**************************************************************************************************************************
// create the sand
void create_dummy_items(ChIrrAppInterface* app, ChSystem& mphysicalSystem, ISceneManager* msceneManager, IVideoDriver* driver)
{
        ChBodySceneNode* mrigidBody;


        int numParticles = 7200;
        int width=63.0 / 2.0;
        int depth = 13;
        int height = numParticles / (width * depth) + 1;
        //int height=11.0 / 2.0;
        //int depth = numParticles / (width * height) + 1;
                //printf("************************* i am adding particle ********************* %d \n", depth);

        int bi = 0;
        for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) {
                        for (int i = 0; i < width; i++) {
                                double Rs = 6.0 / 3.5 * 2;       // scoop 'bounding sphere' rad. for approximate inertia
                                double Ms = 4.0 / 3.0 *CH_C_PI*Rs*Rs*Rs*0.0000026 * massScaleParticles;  // sphere mass (granite)
                                double Js = (2.0/5.0)*Ms*Rs*Rs; // scoop XX YY ZZ inertia, as sphere, just to test things
                                if (bi < numParticles) {
                                        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easySphere(
                                                                                                                                &mphysicalSystem, msceneManager,
                                                                                                                                Ms,
                                                                                                                                ChVector<>(-400 - i * 2.1 * Rs, 20 + j * 2.1f * Rs, -150 + Rs + 2.1f * Rs * k),
                                                                                                                                Rs);
                                        mrigidBody->GetBody()->SetBodyFixed(true);
                                        mrigidBody->GetBody()->SetFriction(0.40);
                                        //mrigidBody->GetBody()->SetRollingFriction( 0.001 );
                                        mrigidBody->GetBody()->SetMass(Ms);
                                        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(Js,Js,Js));
                                        //mrigidBody->addShadowVolumeSceneNode();


                                        video::ITexture* sphereMap = driver->getTexture("../data/rock.jpg");
                                        mrigidBody->setMaterialTexture(0,sphereMap);
                                }
                                bi ++;
                        }
                }
        }
}
//**************************************************************************************************************************
int main(int argc, char* argv[])
{

        // In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
        // global functions are needed.
        ChGlobals* GLOBAL_Vars = DLL_CreateGlobals();

        // Create a ChronoENGINE physical system
        ChSystem mphysicalSystem(60000, 500);
        mphysicalSystem.SetParallelThreadNumber(1);

        // set gravity vector
        mphysicalSystem.Set_G_acc (ChVector<>(0.0, -9800.0, 0.0));

        // Create the Irrlicht visualization (open the Irrlicht device,
        // bind a simple user interface, etc. etc.)

        ChIrrAppInterface application(&mphysicalSystem, L"Excavator",core::dimension2d<u32>(800,600),false);


        /*// This is for GUI tweaking of system parameters..
        MyEventReceiver receiver(&application);
          // note how to add the custom event receiver to the default interface:
        application.SetUserEventReceiver(&receiver);*/

        // Easy shortcuts to add logo, camera, lights and sky in Irrlicht scene:
        ChIrrWizard::add_typical_Logo(application.GetDevice());
        ChIrrWizard::add_typical_Sky(application.GetDevice());
        ChIrrWizard::add_typical_Lights(application.GetDevice(),core::vector3df(-100,370,-700),core::vector3df(-800,800,200),1400.0,1500.0);
        ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-350,270,400),core::vector3df(-300,70,30));

        // choose the solver
        mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR_MULTITHREAD);
        //mphysicalSystem.SetLcpSolverType(ChSystem::LCP_ITERATIVE_GPU);

        // choose the integration type
        //mphysicalSystem.SetIntegrationType(ChSystem::INT_ANITESCU);                                                                                   //?$ zereshk
        mphysicalSystem.SetIntegrationType(ChSystem::INT_TASORA); // for large interpenetrations

        // some solver parameters
        mphysicalSystem.SetIterLCPmaxItersSpeed(50); // constraints can be fulfilled easier the more iterations you use
        mphysicalSystem.SetIterLCPmaxItersStab(2); // only for Tasora integrator
//        ChCollisionModel.SetDefaultSuggestedEnvelope(0.1); //Arman: I don't know why it doesn't work
//        ChCollisionModel.SetDefaultSuggestedMargin  (0.05);
        //mphysicalSystem.SetMaxPenetrationRecoverySpeed(0.2);
        //////////////////////////////////////////////////////////////////////////////////
        // Create the excavator
        video::ITexture* excavatorMap = application.GetVideoDriver()->getTexture("../data/Yellow.png");
        ChBodySceneNode* scoop = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/bucket_mm.obj",
                                                                                0.4,
                                                                                ChVector<>(-447,45,29),
                                                                                QUNIT);

        scoop->GetBody()->SetMass(0.3 * massScale);
        scoop->GetBody()->SetInertiaXX(ChVector<>(290000,360000,260000) * massScale);
        scoop->GetBody()->SetFriction(0.30);
        //scoop ->GetBody()->SetBodyFixed(true);
        //scoop->GetBody()->SetCollide(false);
        scoop->GetBody()->SetCollide(true);
        scoop->GetBody()->GetCollisionModel()->SetFamily(3);
        scoop->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
        scoop->setMaterialTexture(0,excavatorMap);

        ChBodySceneNode* forearm = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/forearm_mm.obj",
                                                                                1.7,
                                                                                ChVector<>(-456,180,29),
                                                                                QUNIT);

        forearm->GetBody()->SetMass(1.7 * massScale);
        forearm->GetBody()->SetInertiaXX(ChVector<>(5000000,800000,5500000) * massScale);
        forearm->GetBody()->SetFriction(0.3);
        //forearm->GetBody()->SetBodyFixed(true);
        forearm->GetBody()->SetCollide(false);
        forearm->GetBody()->GetCollisionModel()->SetFamily(3);
        forearm->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
        forearm->setMaterialTexture(0,excavatorMap);

        ChBodySceneNode* upperarm = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/upperarm_mm.obj",
                                                                                7.2/1.5,
                                                                                ChVector<>(-219,175,29),
                                                                                QUNIT);

        upperarm->GetBody()->SetMass(7.2/1.5 * massScale);
        upperarm->GetBody()->SetInertiaXX(ChVector<>(20000000,60000000,78000000)/1.5 * massScale);
        upperarm->GetBody()->SetFriction(0.3);
        //upperarm->GetBody()->SetBodyFixed(true);
        upperarm->GetBody()->SetCollide(false);
        upperarm->GetBody()->GetCollisionModel()->SetFamily(3);
        upperarm->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
        upperarm->setMaterialTexture(0,excavatorMap);

        ChBodySceneNode* body = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/body_mm.obj",
                                                                                2*17,
                                                                                ChVector<>(-19,80,-5),
                                                                                QUNIT);

        body->GetBody()->SetMass(2*17 * massScale);
        body->GetBody()->SetInertiaXX(2.0 * massScale * ChVector<>(51000000,100000000,71000000));
        body->GetBody()->SetFriction(0.3);
        //body->GetBody()->SetBodyFixed(true);
        body->GetBody()->SetCollide(false);
        body->GetBody()->GetCollisionModel()->SetFamily(3);
        body->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
        body->setMaterialTexture(0,excavatorMap);

        ChBodySceneNode* drive = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/drive_mm.obj",
                                                                                0.2*5,
                                                                                ChVector<>(-20,7,-3),
                                                                                QUNIT);

        drive->GetBody()->SetMass(0.2*5 * massScale);
        drive->GetBody()->SetInertiaXX(0.2*ChVector<>(14000000,29000000,19000000) * massScale);
        drive->GetBody()->SetFriction(0.3);
        //drive->addShadowVolumeSceneNode();
        drive->GetBody()->SetBodyFixed(true);
        drive->GetBody()->SetCollide(false);
        drive->GetBody()->GetCollisionModel()->SetFamily(3);
        drive->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
        drive->setMaterialTexture(0,excavatorMap);

        ChBodySceneNode* tracks = (ChBodySceneNode*)addChBodySceneNode_easyConcaveMesh(
                                                                                application.GetSystem(),
                                                                                application.GetSceneManager(),
                                                                                "../data/tracks_mm.obj",
                                                                                2*5,
                                                                                ChVector<>(-20,7,-3),
                                                                                QUNIT);

        tracks->GetBody()->SetMass(0.2*5 * massScale);
        tracks->GetBody()->SetInertiaXX(0.2*ChVector<>(14000000,29000000,19000000) * massScale);
        tracks->GetBody()->SetFriction(1.0);
        //tracks->addShadowVolumeSceneNode();
        tracks->GetBody()->GetCollisionModel()->BuildModel();
        tracks->GetBody()->SetBodyFixed(true);
        tracks->GetBody()->SetCollide(true);
        tracks->GetBody()->GetCollisionModel()->SetFamily(3);
        tracks->GetBody()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);

        ///////////////////////////////////////////////////////////////////
        video::ITexture* floorMap = application.GetVideoDriver()->getTexture("../data/blu.png");

        // Create the fixed floor
        ChBodySceneNode*        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),                                                                                                application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-300,-36, -235),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(1000,5,410) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(1.0);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);
        //mrigidBody->GetBody()->SetRollingFriction( 0.05 );
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),                                                                                                application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-300,-36, 325),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(1000,5,470) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(1.0);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);
        //mrigidBody->GetBody()->SetRollingFriction( 0.05 );
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),
                                                                                        application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-275,-36, 30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(150,5,120) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(1.0);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);
        //      mrigidBody->GetBody()->SetRollingFriction( 0.05 );
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),
                                                                                        application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(0,-36, 30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(400,5,120) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(1.0);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),
                                                                                        application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-725,-36, 30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(150,5,120) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(1.0);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);
        //mrigidBody->GetBody()->SetRollingFriction( 0.05 );

        // Create a containter for the sand particles
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),
                                                                                        application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-375,-67,30),
                                                                                        ChQuaternion<>(1,0,0,0.0),
                                                                                        ChVector<>(5,90,120) );
        mrigidBody->GetBody()->SetFriction(0.5);
        mrigidBody->GetBody()->SetRot(Q_from_AngAxis(-0.8, VECT_Z));
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);


        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),
                                                                                        application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-600,-66,30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(5,120,120) );
        mrigidBody->GetBody()->SetFriction(0.5);
        mrigidBody->GetBody()->SetRot(Q_from_AngAxis(1.1, VECT_Z));
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);

        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),                                                                                                application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-500,-66,-30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(300,60,5) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(0.5);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);

        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),                                                                                                application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-500,-66, 90),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(300,60,5) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(0.5);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);

        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
                                                                                        application.GetSystem(),                                                                                                application.GetSceneManager(),
                                                                                        10.0,
                                                                                        ChVector<>(-500,-96, 30),
                                                                                        ChQuaternion<>(1,0,0,0),
                                                                                        ChVector<>(300,5,120) );
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->SetFriction(0.5);
        mrigidBody->GetBody()->SetCollide(true);
        mrigidBody->GetBody()->GetCollisionModel()->SetFamily(3);
        mrigidBody->GetBody()->SetInertiaXX(ChVector<>(1000000,1000000,1000000) * massScale);
        mrigidBody->setMaterialTexture(0,floorMap);

        // ------- Create the constraints

        ChSharedPtr<ChLinkEngine> my_motor0(new ChLinkEngine);
        my_motor0->Initialize(tracks->GetBody(), drive->GetBody(), ChCoordsys<>(ChVector<>(-20,7,-3),Q_from_AngAxis(CH_C_PI_2, VECT_Z)) );
        my_motor0->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
        ChFunction_Const * motFunc0 = (ChFunction_Const *)my_motor0->Get_spe_funct();
        motFunc0->Set_yconst(0.0 );
        mphysicalSystem.AddLink(my_motor0);

        ChSharedPtr<ChLinkEngine> my_motor1(new ChLinkEngine);
        my_motor1->Initialize(body->GetBody(), drive->GetBody(), ChCoordsys<>(ChVector<>(-40,50,0),Q_from_AngAxis(CH_C_PI_2, VECT_X)) );
        my_motor1->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
        mphysicalSystem.AddLink(my_motor1);



        ChSharedPtr<ChLinkEngine> my_motor2(new ChLinkEngine);
        my_motor2->Initialize(upperarm->GetBody(), body->GetBody(), ChCoordsys<>(ChVector<>(-90,80,29),Q_from_AngAxis(CH_C_PI_2, VECT_Z)) );
        my_motor2->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
        mphysicalSystem.AddLink(my_motor2);


        ChSharedPtr<ChLinkEngine> my_motor3(new ChLinkEngine);
        my_motor3->Initialize(forearm->GetBody(), upperarm->GetBody(), ChCoordsys<>(ChVector<>(-431,223,29),Q_from_AngAxis(CH_C_PI_2, VECT_Z)) );
        my_motor3->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
        mphysicalSystem.AddLink(my_motor3);


        ChSharedPtr<ChLinkEngine> my_motor4(new ChLinkEngine);
        my_motor4->Initialize(scoop->GetBody(), forearm->GetBody(), ChCoordsys<>(ChVector<>(-485,95,29),Q_from_AngAxis(CH_C_PI_2, VECT_Z)) );
        my_motor4->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
        mphysicalSystem.AddLink(my_motor4);


        // Create some PID objects that will be used to move the bow in 'imposed torque' mode.

        //the contoller for timestep 0.00005s
        /*my_controllerPID1->P=60.0*13000000000*0.5;
        my_controllerPID1->D=60.0*13000000000*0.6*0.12*1.0;
        my_controllerPID1->I=60.0*13000000000*0.5/(0.5*0.5);

        ChSharedPtr<ChControllerPID> my_controllerPID2(new ChControllerPID);
        my_controllerPID2->P=40.0*15000000000*0.5;//0.4 for 0.0025s, 0.2 for 0.005s, 0.5 for 0.001s
        my_controllerPID2->D=40.0*15000000000*0.6*0.12*5.0*0.3;//time constant times 1.5 for 0.0025s, times 2 for 0.005s, times 1 for 0.001s
        my_controllerPID2->I=40.0*15000000000*0.5/(0.5*2.5*0.3);

        ChSharedPtr<ChControllerPID> my_controllerPID3(new ChControllerPID);
        my_controllerPID3->P=80.0*3000000000*0.5;//0.6 for 0.0025s, 0.25 for 0.005s, 0.9 for 0.001s
        my_controllerPID3->D=80.0*3000000000*0.6*0.12*2.4*0.3;//time constant times 1.3 for 0.0025s, times 1.5 for 0.005s, times 1 for 0.001s
        my_controllerPID3->I=80.0*3000000000*0.5/(0.5*1.2*0.3);

        ChSharedPtr<ChControllerPID> my_controllerPID4(new ChControllerPID);
        my_controllerPID4->P=80.0*80000000*0.5;//0.8 for 0.0025s, 0.35 for 0.005s, 0.9 for 0.001s
        my_controllerPID4->D=80.0*80000000*0.6*0.12*2.2*0.3;//time constant times 1.8 for 0.0025s, times 2 for 0.005s, times 1 for 0.001s
        my_controllerPID4->I=80.0*80000000*0.5/(0.5*1.1*0.3);
        */
        ChSharedPtr<ChControllerPID> my_controllerPID1(new ChControllerPID);
        my_controllerPID1->P=1.5*13000000000*0.5 * massScale;
        my_controllerPID1->D=1.5*13000000000*0.6*0.12*5.0*7.0 * massScale;
        my_controllerPID1->I=1.5*13000000000*0.5/(0.5*2.5*7.0) * massScale;

        ChSharedPtr<ChControllerPID> my_controllerPID2(new ChControllerPID);
        my_controllerPID2->P=0.4*15000000000*0.5 * massScale;//0.4 for 0.0025s, 0.2 for 0.005s, 0.5 for 0.001s
        my_controllerPID2->D=0.4*15000000000*0.6*0.12*5.0*1.0*7.0 * massScale;//time constant times 1.5 for 0.0025s, times 2 for 0.005s, times 1 for 0.001s
        my_controllerPID2->I=0.4*15000000000*0.5/(0.5*2.5*1.0*7.0) * massScale;

        ChSharedPtr<ChControllerPID> my_controllerPID3(new ChControllerPID);
        my_controllerPID3->P=0.9*3000000000*0.5 * massScale;//0.6 for 0.0025s, 0.25 for 0.005s, 0.9 for 0.001s
        my_controllerPID3->D=0.9*3000000000*0.6*0.12*2.4*1.0*5.0 * massScale;//time constant times 1.3 for 0.0025s, times 1.5 for 0.005s, times 1 for 0.001s
        my_controllerPID3->I=0.9*3000000000*0.5/(0.5*1.2*1.0*5.0) * massScale;

        ChSharedPtr<ChControllerPID> my_controllerPID4(new ChControllerPID);
        my_controllerPID4->P=0.9*80000000*0.5 * massScale;//0.8 for 0.0025s, 0.35 for 0.005s, 0.9 for 0.001s
        my_controllerPID4->D=0.9*80000000*0.6*0.12*2.2*1.0*5.0 * massScale;//time constant times 1.8 for 0.0025s, times 2 for 0.005s, times 1 for 0.001s
        my_controllerPID4->I=0.9*80000000*0.5/(0.5*1.1*1.0*5.0) * massScale;

        //////////////////////////////////////////////////////////////////////////////////
        // THE TIME INTEGRATION CYCLE
        //

        int padnumber=0;
        int save=0;
        double dT = .0001;//.00001;//.0005;
        double periodExtension = .0005 / dT;
    int savestep = 10 * periodExtension;
        // initialize the memory of contoller positions
        double mem_control[]={0,0,0,0};

        // relative rotations for the motion controller
        // rot_motor1, rot_motor2, rot_motor3, rot_motor4, tolerance, numsteps, ....
        double rotation_control[]=      {0.0,-0.25,0.0,0.0,0.01,800/2.5 * periodExtension / 3.0,//move up
                                0.0,0.0,0.0,0.0,0.01,1500/2.5 * periodExtension,//waiting
                                0.0,0.45,-0.9,-0.9,0.01,2000/2.5 * periodExtension,
                                0.0,0.03,1.1,1.3,0.01,6000/2.5 * periodExtension,//digging
                                0.0,-0.06,0.0,0.0,0.01,250/2.5 * periodExtension,//move up//
                                0.6,-0.21,0.0,0.2,0.01,2000/2.5 * periodExtension,//turn
                                0.0,-0.24,-0.8,-1.3,0.01,4000/2.5 * periodExtension,//release
                                -0.6,0.47,-0.3,-0.2,0.01,2000/2.5 * periodExtension,//turn back
                                0.0,0.05,1.1,1.25,0.01,6000/2.5 * periodExtension,//digging
                                0.0,-0.11,0.0,0.0,0.01,400/2.5 * periodExtension,//move up
                                0.6,-0.15,0.0,0.2,0.01,2000/2.5 * periodExtension,//turn
                                0.0,-0.27,-0.8,-1.3,0.01,4000/2.5 * periodExtension,//release
                                -0.6,0.50,-0.3,-0.2,0.01,2000/2.5 * periodExtension,//turn back
                                0.0,0.05,1.1,1.25,0.01,6000/2.5 * periodExtension,//digging
                                0.0,-0.11,0.0,0.0,0.01,400/2.5 * periodExtension,//move up
                                0.6,-0.15,0.0,0.2,0.01,2000/2.5 * periodExtension,//turn
                                0.0,-0.24,-0.8,-1.3,0.01,4000/2.5 * periodExtension//release
                                //0.0,-0.25,0.0,-0.5,0.01,1200/2.5,
                                //0.0,-0.05,0.8,-0.4,0.01,120/2.5,//throughing
                                //0.0,0.0,0.0,0.0,0.01,10000/2.5//waiting
                                //0.0,0.0,-0.3,0.6,0.01,3000/5,
                                //0.0,0.8,-0.2,-0.4,0.01,3000/5,//push up
                                //0.0,-0.3,0.0,0.0,0.01,1000/5
                                };

        // initialize variables for motioncontrol
        int myswitch=0;
        int stage=0;
        int step=0;
        int num_step=1;
        int num_stages=0;
        // create variables for the file saving step
        ChQuaternion<> bodyRot;
        ChVector<> bodyAngs;

        int cumulStep = 0;
        int startWrite = rotation_control[5] + rotation_control[11];
        int endWrite = startWrite + rotation_control[17] + rotation_control[23] + rotation_control[29];

        clock_t startTime, endTime;
        startTime = clock();

        time_t startTimeBeta = time(NULL);

        fstream outFileTime;
        outFileTime.open("simTimeDVI_CPU.txt", ios::out);

        ///dummy items
        create_dummy_items(&application, mphysicalSystem, application.GetSceneManager(), application.GetVideoDriver());


        while(application.GetDevice()->run() && cumulStep < endWrite) {
                step++;
                cumulStep++;
                //printf("step:%d/%d ", cumulStep, endWrite);

                if ((cumulStep - int(rotation_control[5])) % (200 * int(periodExtension)) == 1 && cumulStep > rotation_control[5] && cumulStep < startWrite) {
                        //printf(" hey you: %d %d %d %d \n", cumulStep, int(rotation_control[5]) % 500, rotation_control[5], startWrite);
                        // Create some more particles
                        create_some_falling_items(&application, mphysicalSystem, application.GetSceneManager(), application.GetVideoDriver());
                }

                // create loop motion
                num_stages=(int)(sizeof(rotation_control)/(6*sizeof(rotation_control[0])));
                stage=(stage+myswitch)%(num_stages+1);
                // control motion
                myswitch=MotionControl(rotation_control, stage, num_stages, step, mem_control, my_motor1, my_motor2, my_motor3, my_motor4);
                if (myswitch==1) step=0;

                // update the PID objects given the setpoint error = user_desired_rotation - actual_rotation
                // and get the needed torque:
                double torque1 = my_controllerPID1->Get_Out( STATIC_motor1 - my_motor1->Get_mot_rot(), mphysicalSystem.GetChTime());
                //torque1 = fast_sign(torque1)*MIN(3000000000000,fabs(torque1));
                ChFunction_Const * torqFunc1 = (ChFunction_Const*)my_motor1->Get_tor_funct();
                torqFunc1->Set_yconst(torque1);
                //ChFunction_Const * torqFunc1 = (ChFunction_Const*)my_motor1->Get_rot_funct();
                //torqFunc1->Set_yconst(STATIC_motor1);

                double torque2 = my_controllerPID2->Get_Out( STATIC_motor2 - my_motor2->Get_mot_rot(), mphysicalSystem.GetChTime());
                //torque2 = fast_sign(torque2)*MIN(3000000000000,fabs(torque2));
                ChFunction_Const * torqFunc2 = (ChFunction_Const*)my_motor2->Get_tor_funct();
                torqFunc2->Set_yconst(torque2);
                //ChFunction_Const * torqFunc2 = (ChFunction_Const*)my_motor2->Get_rot_funct();
                //torqFunc2->Set_yconst(STATIC_motor2);

                double torque3 = my_controllerPID3->Get_Out( STATIC_motor3 - my_motor3->Get_mot_rot(), mphysicalSystem.GetChTime());
                //torque3 = fast_sign(torque3)*MIN(2000000000000,fabs(torque3));
                ChFunction_Const * torqFunc3 = (ChFunction_Const*)my_motor3->Get_tor_funct();
                torqFunc3->Set_yconst(torque3);
                //ChFunction_Const * torqFunc3 = (ChFunction_Const*)my_motor3->Get_rot_funct();
                //torqFunc3->Set_yconst(STATIC_motor3);

                double torque4 = my_controllerPID4->Get_Out( STATIC_motor4 - my_motor4->Get_mot_rot(), mphysicalSystem.GetChTime());
                //torque4 = fast_sign(torque4)*MIN(1000000000000,fabs(torque4));
                ChFunction_Const * torqFunc4 = (ChFunction_Const*)my_motor4->Get_tor_funct();
                torqFunc4->Set_yconst(torque4);
                //ChFunction_Const * torqFunc4 = (ChFunction_Const*)my_motor4->Get_rot_funct();
                //torqFunc4->Set_yconst(STATIC_motor4);

                if ( (cumulStep == 1) || (cumulStep % (10 * int(periodExtension)) == 0 && cumulStep >= startWrite && cumulStep < endWrite) ) {
                        PrintToFileForces(my_motor1, torque1, "my_motor1.txt", cumulStep, dT, startWrite, endWrite);
                        PrintToFileForces(my_motor2, torque2, "my_motor2.txt", cumulStep, dT, startWrite, endWrite);
                        PrintToFileForces(my_motor3, torque3, "my_motor3.txt", cumulStep, dT, startWrite, endWrite);
                        PrintToFileForces(my_motor4, torque4, "my_motor4.txt", cumulStep, dT, startWrite, endWrite);

                        PrintToFile_BodyForces(scoop->GetBody(), "scoopForces.txt", cumulStep);



                        //PrintToFileForces(my_motor1, torque1, "C:/Users/Arman/Documents/cat_proj/Data/my_motor1.txt", cumulStep, dT, startWrite, endWrite);
                        //PrintToFileForces(my_motor2, torque2, "C:/Users/Arman/Documents/cat_proj/Data/my_motor2.txt", cumulStep, dT, startWrite, endWrite);
                        //PrintToFileForces(my_motor3, torque3, "C:/Users/Arman/Documents/cat_proj/Data/my_motor3.txt", cumulStep, dT, startWrite, endWrite);
                        //PrintToFileForces(my_motor4, torque4, "C:/Users/Arman/Documents/cat_proj/Data/my_motor4.txt", cumulStep, dT, startWrite, endWrite);

                }


                // do time integration
                mphysicalSystem.DoStepDynamics( dT );
                //time
                if (cumulStep % savestep == 1) {outFileTime<<"realTime:"<<dT * cumulStep<<", simulation time by far (s): "<<time(NULL) - startTimeBeta<<endl;}

                //irrlicht graphic  window, true: enable, false: disable
                if (false) {
                        //**************************** create the irrlich visualization
                        if (cumulStep % savestep == 1) {
                                application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));
                                application.DrawAll();
                        }
                        // save screenshots to file, to make an AVI later

                        if (cumulStep % savestep == 1) {
                                padnumber++;
                        video::IImage* image = application.GetDevice()->getVideoDriver()->createScreenShot();
                        char filename1[100];
                                sprintf(filename1, "screenshots3/screenshot%04d.png", padnumber);

                        GetLog() << "\n saving screenshot: " << filename1;

                                if (image)
                                application.GetDevice()->getVideoDriver()->writeImageToFile(image, filename1);
                        image->drop();
                        }
                        // write a file with all body positions for the use in rendering programs
                        char pad_number[100];
                        sprintf(pad_number, "%d", (padnumber+1000000));

                        char filename[100];
                        if (padnumber<200000) {
                                sprintf(filename, "excavatorSim_data/Pos%s.dat", pad_number+1);
                        }
                        else if (padnumber<400000) {
                                sprintf(filename, "excavatorSim_data2/Pos%s.dat", pad_number+1);
                        }
                        else {
                                sprintf(filename, "excavatorSim_data3/Pos%s.dat", pad_number+1);
                        }
                        ///ChStreamOutAsciiFile data_spheres_positions(filename);
                        //std::vector<ChBody*>::iterator abody = mphysicalSystem.Get_bodylist()->begin();
                        //while (abody != mphysicalSystem.Get_bodylist()->end())
                        //      {
                        //              ChBody* bpointer = (*abody);
                        //              bodyRot = bpointer->GetRot();
                        //              bodyAngs = bodyRot.Q_to_NasaAngles();

                        //              data_spheres_positions << bpointer->GetPos().x << ", ";
                        //              data_spheres_positions << bpointer->GetPos().y << ", ";
                        //              data_spheres_positions << bpointer->GetPos().z << ", ";
                        //              data_spheres_positions << bodyAngs(0) << ", ";
                        //              data_spheres_positions << bodyAngs(1) << ", ";
                        //              data_spheres_positions << bodyAngs(2) << ", ";
                        //              data_spheres_positions << "\n";
                        //              abody++;
                        //      }
                        if (cumulStep % savestep == 1)
                                application.GetVideoDriver()->endScene();
                }
        }
        //if (application.GetDevice()->run()) application.GetVideoDriver()->endScene();
        endTime = clock();
        outFileTime.close();
        double simulationTime = (endTime - startTime) / double(CLOCKS_PER_SEC);
        printf("simulation time (ms) %f\n", 1000 * simulationTime);

        // Remember this at the end of the program, if you started
        // with DLL_CreateGlobals();
        DLL_DeleteGlobals();
        return 0;


}

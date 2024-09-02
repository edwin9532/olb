#include "olb3D.h"
#include "olb3D.hh"
#include "cmath"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q7<>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

// Parameters for the simulation setup
const T lx0   = 10.0;    // length of channel
const T ly0   = 10.0;     // height of channel
const T lz0   = 10.0;     // width of channel
//const int N = 30;         // resolution of the model
//const int M = 30;         // resolution of the model (time)
const T maxPhysT = 6.;  // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  //todo el dominio pasa a ser de material 2
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );

    // Agregar una celda central con un nuevo material
  Vector<T, 3> center((lx0/2), (ly0/2), (lz0/2));  // Posición central
  Vector<T, 3> cellSize(converter.getConversionFactorLength(),
                        converter.getConversionFactorLength(),
                        converter.getConversionFactorLength());  // Tamaño de una celda

  // Crear un cuboide que representa la celda central
  IndicatorCuboid3D<T> centralCell(2.*cellSize, center);

  // Renombrar la celda central como material 3
  superGeometry.rename(1, 3, centralCell);
 
  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T,DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  

  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3});

  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  // Material=2 -->bounce back
  setBounceBackBoundary(sLattice, superGeometry, 2);

  // Setting of the boundary conditions


  // Initial conditions
  AnalyticalConst3D<T,T> ux( 0. );
  AnalyticalConst3D<T,T> uy( 0. );
  AnalyticalConst3D<T,T> uz( 0. );

  AnalyticalConst3D<T,T> rho( 1. );

  AnalyticalComposed3D<T,T> u( ux,uy,uz );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rho, u );
  sLattice.iniEquilibrium( bulkIndicator, rho, u );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(UnitConverter<T, DESCRIPTOR> const& converter,
                       SuperLattice<T, DESCRIPTOR>& sLattice, unsigned int iT,
                       SuperGeometry<T, 3>& superGeometry)
{
    OstreamManager clout(std::cout, "setBoundaryValues");
    auto changeTime = converter.getLatticeTime(2.);
    // Definir el tiempo en el que quieres cambiar la densidad
    // por ejemplo, 20 unidades de tiempo

    // Cambiar la densidad en el material 3 cuando se alcance el tiempo especificado
    if (iT < changeTime) {
        sLattice.setProcessingContext(ProcessingContext::Evaluation);
        T omega = 0.4;
        T newDensity = 0.3*sin(omega*iT)+1; // Ajusta este valor según sea necesario
        AnalyticalConst3D<T, T> rho(newDensity);
        sLattice.defineRho(superGeometry, 3, rho);

        //clout << "Densidad cambiada a " << newDensity << " en el paso de tiempo " << iT << std::endl;

        sLattice.setProcessingContext(ProcessingContext::Simulation);
    }
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter,
                 BlockReduction3D2D<T>& planeReduction,
                 int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "bstep3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int  vtkIter  = converter.getLatticeTime( 0.1 );
  const int  statIter = converter.getLatticeTime( 0.1 );
  //const int  saveIter = converter.getLatticeTime( 1. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    vtmWriter.write( iT );
    planeReduction.update();
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = converter.getCharPhysVelocity() * 3./2.;
    jpeg_Param.minValue = 0.0;
    jpeg_Param.fullScreenPlot = true;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Saves lattice data
  //if ( iT%( saveIter/2 )==0 && iT>0 ) {
  //  clout << "Checkpointing the system at t=" << iT << std::endl;
  //  sLattice.save( "bstep3d.checkpoint" );
  //  // The data can be reloaded using
  //  //     sLattice.load("bstep3d.checkpoint");
  //}
}

int main( int argc, char* argv[] )
{
  int N = std::atoi(argv[1]);         // resolution of the model         :30
  int M = std::atoi(argv[2]);         // resolution of the model (time)  :30
  
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverter<T,DESCRIPTOR> converter(
    (T)   1./N,     // physDeltaX: spacing between two lattice cells in __m__
    (T)   1./(M*N), // physDeltaT: time step in __s__
    (T)   1.,       // charPhysLength: reference length of simulation geometry
    (T)   1.,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./200.,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.        // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep3d");

  // === 2nd Step: Prepare Geometry ===
  STLreader<T> stlReader( "10000mm-XY.stl", converter.getConversionFactorLength(), 0.001 );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  // Set up persistent measuring functors for result extraction
  SuperLatticePhysPressure3D<T, DESCRIPTOR> velocity( sLattice, converter );
  //SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperEuklidNorm3D<T> normVel( velocity );

  BlockReduction3D2D<T> planeReduction(
    normVel,
    Hyperplane3D<T>().centeredIn(cuboidGeometry.getMotherCuboid()).normalTo({0,0,1}),
    600,
    BlockDataSyncMode::ReduceOnly);
  
  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    //getResults( sLattice, converter, planeReduction, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printShortSummary();
}


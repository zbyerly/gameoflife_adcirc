/**
 * We need to include typemaps first to avoid problems with Intel
 * MPI's C++ bindings (which may collide with stdio.h's SEEK_SET,
 * SEEK_CUR etc.).
 */
#include <libgeodecomp/communication/typemaps.h>
#include <libgeodecomp/communication/mpilayer.h>
#include <libgeodecomp/parallelization/serialsimulator.h>
#include <libgeodecomp/parallelization/stripingsimulator.h>

// Added by ZDB
#include <libgeodecomp/storage/multicontainercell.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <math.h>
// ------------

#include <boost/assign/std/vector.hpp>

#include <libgeodecomp/config.h>
#include <libgeodecomp/geometry/stencils.h>
#include <libgeodecomp/io/bovwriter.h>
#include <libgeodecomp/io/ppmwriter.h>
#include <libgeodecomp/io/simplecellplotter.h>
#include <libgeodecomp/io/simpleinitializer.h>
#include <libgeodecomp/io/tracingwriter.h>
#include <libgeodecomp/loadbalancer/oozebalancer.h>
#include <libgeodecomp/loadbalancer/tracingbalancer.h>
#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/storage/image.h>

using namespace boost::assign;
using namespace LibGeoDecomp;

// Added by ZDB
FloatCoord<2> origin;
FloatCoord<2> quadrantDim;
// ------------

// Each instance represents one subdomain of an ADCIRC unstructured grid
class DomainCell
{
public:
    class API: 
        public::APITraits::HasCustomRegularGrid,
        public::APITraits::HasUnstructuredGrid,
        public::APITraits::HasPointMesh
    {
    public:
        inline FloatCoord<2> getRegularGridSpacing()
        {
            return quadrantDim;
        }
        
        inline FloatCoord<2> getRegularGridOrigin()
        {
            return origin;
        }
    };

    DomainCell(const LibGeoDecomp::FloatCoord<2>& center = FloatCoord<2>(), int id = 0, int wetdry = 0, int numResidentNodes = 0) :
        center(center),
        id(id),
        wetdry(wetdry),
        numResidentNodes(numResidentNodes)        
    {}

    template<typename NEIGHBORHOOD>
    void update(const NEIGHBORHOOD& hood, int nanoStep);

    void pushNeighborNode(const int neighborID)
    {
        if (std::count(neighboringNodes.begin(), neighboringNodes.end(), neighborID) == 0) {
            neighboringNodes << neighborID;
        }
    }

    void pushResidentNode(const int resNodeID)
    {
        this->residentNodes[this->numResidentNodes] = resNodeID;
        this->numResidentNodes++;        
    }
    

    const LibGeoDecomp::FloatCoord<2>& getPoint() const
    {
        return center;
    }

    // fixme: I don't like that I have to specify this despite not needing it
    std::vector<LibGeoDecomp::FloatCoord<2> > getShape() const
    {
        std::vector<LibGeoDecomp::FloatCoord<2> > ret;
        ret << center + FloatCoord<2>( 0.000, -0.001);
        ret << center + FloatCoord<2>(-0.001,  0.000);
        ret << center + FloatCoord<2>( 0.000,  0.001);
        ret << center + FloatCoord<2>( 0.001,  0.000);

        return ret;
    }

    LibGeoDecomp::FloatCoord<2> center;
    int id;
    int wetdry;
  // TODO: change these fixedarrays to something else
  std::vector<int> neighboringNodes;
  std::vector<int> residentNodes;
  std::vector<int> residentNodesAlive;
  int numResidentNodes;
    
};
// ContainerCell translates between the unstructured grid and the
// regular grid currently required by LGD
typedef ContainerCell<DomainCell, 1000> ContainerCellType;

// type definition required for callback functions below
typedef LibGeoDecomp::TopologiesHelpers::Topology<2, false, false, false > TopologyType;
typedef LibGeoDecomp::Grid<ContainerCellType, TopologyType> GridType;
typedef LibGeoDecomp::CoordMap<ContainerCellType, GridType> BaseNeighborhood;
typedef LibGeoDecomp::NeighborhoodAdapter<BaseNeighborhood, 2, DomainCell> Neighborhood;

const Neighborhood *neighborhood;
DomainCell *domainCell;

template<typename NEIGHBORHOOD>
void DomainCell::update(const NEIGHBORHOOD& hood, int nanoStep)
{
    std::cerr << "Entering c subroutine update\n";
    // init callback, by which the Fortran code can retrieve data
    domainCell = this;
    neighborhood = &hood;
    int numNodeNeighbors = neighboringNodes.size();

//    update_node_(&id, &wetdry, &neighboringNodes[0], &numNodeNeighbors);
    //outputBorderNodes is a vector of ints which are the "outgoing"
    //border nodes, in *local* indices
    int numResidentNodes=this->numResidentNodes;
}



class ADCIRCInitializer : public SimpleInitializer<ContainerCellType>
{
public:
    //typedef typename ContainerCellType::Cargo Cargo;
    typedef GridBase<ContainerCellType, 2> GridType;
    using SimpleInitializer<ContainerCellType>::dimensions;    

    ADCIRCInitializer(const std::string& meshDir, const int steps) :
        SimpleInitializer<ContainerCellType>(Coord<2>(), steps),
        meshDir(meshDir)
    {
        std::cerr << "determining grid dimensions\n";
        
        determineGridDimensions();
    }
    
    virtual void grid(GridType *grid)
    {
        std::cerr << "Entering c subroutine grid\n";
        std::ifstream fort80File;
        
        int numberOfDomains;
        int numberOfElements;
        int numberOfPoints;
        std::vector<ownerTableEntry> ownerTable;
        std::vector<neighborTable> myNeighborTables;
        
        openfort80File(fort80File);
        readfort80(fort80File, &numberOfDomains, &numberOfElements, &numberOfPoints, &ownerTable);


        std::vector<std::vector<int> > neighboringDomains;        
        std::vector<FloatCoord<2> > centers;


        // clear grid:
        CoordBox<2> box = grid->boundingBox();
        for (CoordBox<2>::Iterator i = box.begin(); i != box.end(); ++i) {
            grid->set(*i, ContainerCellType());
        }


        // piece together domain node cells:
        for (int i=0; i< numberOfDomains; i++){            
            neighborTable myNeighborTable;
            std::ifstream fort18File;
            int numberOfNeighbors;
            openfort18File(fort18File, i);
            readfort18(fort18File, &numberOfNeighbors, &i, &myNeighborTable);            
//            std::cerr << "Domain: " << i << " number of neighbors: " << numberOfNeighbors << "\n";
            myNeighborTables.push_back(myNeighborTable);

            //Read fort.14 file for each domain
            int numberOfPoints;
            int numberOfElements;
            std::ifstream fort14File;
            openfort14File(fort14File, i);
            readFort14Header(fort14File, &numberOfElements, &numberOfPoints);
            std::vector<FloatCoord<3> > points;
            readFort14Points(fort14File, &points, numberOfPoints);            
            FloatCoord<2> center = determineCenter(&points);
            center /= numberOfPoints;
            
            centers.push_back(center);

            // Load neighboringDomains with myNeighborTables
            std::vector<int> neighbors;
            for (int j=0; j<numberOfNeighbors; j++){
                neighbors.push_back(myNeighborTable.myNeighbors[j].neighborID);
            }
            // I'm not sure if I need to do this.
            neighboringDomains.push_back(neighbors);

            int nodeID = i;
            DomainCell node(center, nodeID);
            
            for (int j=0; j<numberOfNeighbors; j++)
            {
                node.pushNeighborNode(neighbors[j]);
            }
            
            for (int j=0; j<ownerTable.size(); j++)
            {
                if (ownerTable[j].ownerID == nodeID)
                {
                    std::cerr << "pushing to owner " << nodeID  << "\n";
                    node.pushResidentNode(ownerTable[j].globalID);
                }
            }
            
                
            
            FloatCoord<2> gridCoordFloat = (node.center - minCoord) / quadrantDim;
            Coord<2> gridCoord(gridCoordFloat[0], gridCoordFloat[1]);

            ContainerCellType container = grid->get(gridCoord);
            container.insert(node.id, node);
            grid->set(gridCoord, container);            
        }        
    }
    
    
    
    
    
private:
    std::string meshDir;
    double maxDiameter;
    FloatCoord<2> minCoord;
    FloatCoord<2> maxCoord;
    
    struct neighbor
    {
        int neighborID;
        std::vector<int> sendNodes;
        std::vector<int> recvNodes;
    };

    struct neighborTable 
    {
        std::vector<neighbor> myNeighbors;
    };

    struct ownerTableEntry
    {
        int globalID;
        int localID;
        int ownerID;
    };
    

    void determineGridDimensions()
    {
        std::cerr << "Entering c subroutine determineGridDimensions\n";
        std::ifstream fort80File;
        
        int numberOfDomains;
        int numberOfElements;
        int numberOfPoints;
        std::vector<ownerTableEntry> ownerTable;
        std::vector<neighborTable> myNeighborTables;
        openfort80File(fort80File);
        readfort80(fort80File, &numberOfDomains, &numberOfElements, &numberOfPoints, &ownerTable);
        
        std::vector<FloatCoord<2> > centers;

        for (int i=0; i< numberOfDomains; i++){
            neighborTable myNeighborTable;
            std::ifstream fort18File;
            int numberOfNeighbors;
            openfort18File(fort18File, i);
            readfort18(fort18File, &numberOfNeighbors, &i, &myNeighborTable);            
            myNeighborTables.push_back(myNeighborTable);

            //Read fort.14 file for each domain
            int numberOfPoints;
            int numberOfElements;
            std::ifstream fort14File;
            openfort14File(fort14File, i);
            readFort14Header(fort14File, &numberOfElements, &numberOfPoints);
            std::vector<FloatCoord<3> > points;
            readFort14Points(fort14File, &points, numberOfPoints);            
            FloatCoord<2> center = determineCenter(&points);
            center /= numberOfPoints;
            
            centers.push_back(center);
        }        

        maxDiameter = determineMaximumDiameter(&centers, myNeighborTables);
        minCoord = FloatCoord<2>(
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max());
        maxCoord = FloatCoord<2>(
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max());

        for (int i = 0; i < numberOfDomains; ++i) {
            FloatCoord<2> p(centers[i][0], centers[i][1]);

            minCoord = minCoord.min(p);
            maxCoord = maxCoord.max(p);
        }

        origin = minCoord;
        // add a safety factor for the cell spacing so we can be sure
        // neighboring elements are never more than 1 cell apart in the grid:
        quadrantDim = FloatCoord<2>(maxDiameter * 2.0, maxDiameter * 2.0);
        FloatCoord<2> floatDimensions = (maxCoord - minCoord) / quadrantDim;
        dimensions = Coord<2>(
            ceil(floatDimensions[0]),
            ceil(floatDimensions[1]));

        
        std::cout << "geometry summary:\n"
                  << "  minCoord: "    << minCoord    << "\n"
                  << "  maxCoord: "    << maxCoord    << "\n"
                  << "  maxDiameter: " << maxDiameter << "\n"
                  << "  dimensions: "  << dimensions  << "\n"
                  << "\n";
    }
    

    FloatCoord<2> determineCenter(std::vector<FloatCoord<3> > *points) 
    {
        FloatCoord<2> center;
//        std::cerr << "points->size() = " << points->size() << "\n";
        for (int i=0; i < points->size(); i++){
//            std::cerr << "points[0][" << i << "][0] = " << points[0][i][0] << "\n";
//            std::cerr << "points[1][" << i << "] = " << points[1][i] << "\n";
//            std::cerr << "points[" << i << "][0] = " << points[i][0] << "\n";
//            std::cerr << "points[" << i << "][1] = " << points[i][1] << "\n";
            FloatCoord<2> coord(points[0][i][0],points[0][i][1]);
            center += coord;
        }
        return center;
    }
    
    double determineMaximumDiameter(
        const std::vector<FloatCoord<2> >* points,
        const std::vector<neighborTable> myNeighborTables)
    {
        double maxDiameter = 0;

        int numPoints = points->size();
//        std::cerr << "numPoints = " << numPoints << "\n";
        for (int point = 0; point < numPoints; ++point) {
            int numNeighbors = myNeighborTables[point].myNeighbors.size();
//            std::cerr << "domain: " <<point << " numNeighbors = " << numNeighbors << "\n";            
            for (int i=0; i<numNeighbors; i++){
                int neighborID = myNeighborTables[point].myNeighbors[i].neighborID;
                double dist = getDistance(points[0][point],points[0][neighborID]);
//                std::cerr << "neighborID = " << neighborID << " distance = " << dist << "\n";
                maxDiameter = std::max(maxDiameter,dist);
            }
        }
        
        return maxDiameter;
    }

    double getDistance(FloatCoord<2> a, FloatCoord<2> b)
    {
        FloatCoord<2> c = a-b;
        return sqrt(c[0]*c[0]+c[1]*c[1]);
    }
    
        
    void openfort80File(std::ifstream& meshFile)
    {
        std::cerr << "Entering c subroutine openfort80File\n";
        std::string meshFileName = meshDir;
        meshFileName.append("/fort.80");
        meshFile.open(meshFileName.c_str());
        if (!meshFile) {
            throw std::runtime_error("could not open fort.80 file "+meshFileName);
        } else {
            std::cerr << "opened fort.80 file " << meshFileName << " successfully.\n";
        }
        
    }
        
    void readfort80(std::ifstream& meshFile, int *numberOfDomains, int *numberOfElements, int *numberOfPoints, std::vector<ownerTableEntry> *ownerTable)
    {
        std::cerr << "Entering c subroutine readfort80\n";
        std::string buffer(1024, ' ');
        //Discard first 3 lines:
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
            
        meshFile >> *numberOfElements;
        meshFile >> *numberOfPoints;
            
        //Discard rest of line
        std::getline(meshFile, buffer);
            
        meshFile >> *numberOfDomains;
            
//        std::cerr << "number of elements: " << *numberOfElements << "\n"
//                  << "number of points: " << *numberOfPoints << "\n"
//                  << "number of domains: " << *numberOfDomains << "\n";    
            
        //Discard rest of the line:
        std::getline(meshFile, buffer);

        //Discard next 8 lines:
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        std::getline(meshFile, buffer);
        
//        std::cerr << buffer << "\n";

        for (int domain = 0; domain < *numberOfDomains; ++domain) {
            int buf;
            int numNodes;
            meshFile >> buf;
            if (buf != domain) {
                throw std::runtime_error("buf does not match domain");
            }
            meshFile >> numNodes;
            //Discard rest of the line
            std::getline(meshFile, buffer);
            
//            std::cerr << "Domain: " << domain << " numNodes: " << numNodes << "\n";
            
            for (int node = 0; node < numNodes; ++node) {
                int nodeNum;
                meshFile >> nodeNum;
                
//                std::cerr << "Domain: " << domain << " Node: " << nodeNum << "\n";
            }
            
            // throw away the rest of the line
            std::getline(meshFile, buffer);
        }
        // Throw away another line
        std::getline(meshFile, buffer);
        for (int node = 0; node < *numberOfPoints; node++) {
            ownerTableEntry thisOwnerTableEntry;
            int global_label;
            int local_label;
            int owner;
            meshFile >> global_label;
            meshFile >> owner;
            meshFile >> local_label;
            
            thisOwnerTableEntry.globalID=global_label;
            thisOwnerTableEntry.ownerID=owner;
            thisOwnerTableEntry.localID=local_label;
            
            std::getline(meshFile, buffer);//discard rest of line
            
            ownerTable->push_back(thisOwnerTableEntry);
        }            
    }    

    void openfort18File(std::ifstream& meshFile, int domainID)
    {
        std::cerr << "Entering c subroutine openfort18File\n";
        std::string meshFileName = meshDir;
        meshFileName.append("/PE");
        std::stringstream buf;
        buf.width(4);
        buf.fill('0');
        buf << domainID;
        meshFileName.append(buf.str());
        meshFileName.append("/fort.18");
//        std::cerr << meshFileName << "\n";
        meshFile.open(meshFileName.c_str());
        if (!meshFile) {
            throw std::runtime_error("could not open fort.18 file"+meshFileName);
        }    
    }


    void readfort18(std::ifstream& meshFile, int *numberOfNeighbors, const int *domainID, neighborTable *myNeighborTable)
    {
        std::cerr << "Entering c subroutine readfort18\n";
        int numberOfElements;
        int numberOfNodes;
        int domainIDFromFile;
        int numberOfResNodes;
        std::vector<int> residentNodes;        

        std::string buffer(1024, ' ');
        //Discard first line:
        std::getline(meshFile, buffer);

        //Discard first three entries of line 2:
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> buffer;
        
        meshFile >> numberOfElements;
        std::getline(meshFile, buffer); //Discard the rest of the line

        for (int i=0; i<numberOfElements; i++){std::getline(meshFile, buffer);}

        //Discard first three entries of next line:
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> buffer;        
        meshFile >> numberOfNodes;
        std::getline(meshFile, buffer); //Discard the rest of the line
        for (int i=0; i<numberOfNodes; i++){std::getline(meshFile, buffer);}

        //NFLUXF line: manually discard for now.  May need to do
        //something more robust later.
        std::getline(meshFile, buffer);

        //Discard first three entries of next line:        
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> buffer;        
        int NETA;
        meshFile >> NETA;
        std::getline(meshFile, buffer); //Discard the rest of the line
        for (int i=0; i<NETA; i++){std::getline(meshFile, buffer);}

        //Discard first three entries of next line:        
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> buffer;        
        int NSTAE;
        meshFile >> NSTAE;
        std::getline(meshFile, buffer); //Discard the rest of the line
        for (int i=0; i<NSTAE; i++){std::getline(meshFile, buffer);}

        //Discard first three entries of next line:        
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> buffer;        
        int NSTAV;
        meshFile >> NSTAV;
        std::getline(meshFile, buffer); //Discard the rest of the line
        for (int i=0; i<NSTAV; i++){std::getline(meshFile, buffer);}
        
        //Discard NSTAM and NSTAC lines manually
        std::getline(meshFile, buffer);        
        std::getline(meshFile, buffer);
        
        //RES NODE line:
        meshFile >> buffer;
        meshFile >> buffer;
        meshFile >> domainIDFromFile;
        meshFile >> numberOfResNodes;
        
        for (int i=0; i<numberOfResNodes; i++)
        {
            int node;
            
            meshFile >> node;
            residentNodes.push_back(node);
        }        
        std::getline(meshFile, buffer); //Discard remainder of the
                                        //line.

        //COMM PE line
        meshFile >> buffer; //COMM
        meshFile >> buffer; //PE
        meshFile >> *numberOfNeighbors;
        std::getline(meshFile, buffer);//discard rest of the line

        for (int i=0; i<*numberOfNeighbors; i++){
            neighbor neighbor;
            int numberOfRecvNodes;
            meshFile >> buffer; //RECV
            meshFile >> buffer; //PE
            meshFile >> neighbor.neighborID;
            meshFile >> numberOfRecvNodes;
            std::getline(meshFile, buffer);//discard rest of the line

            //Assemble arrays of nodes to be received
            for (int j=0; j<numberOfRecvNodes; j++){
                int receiveNode;
                meshFile >> receiveNode;   
                neighbor.recvNodes.push_back(receiveNode);
            }
            std::getline(meshFile, buffer);//discard rest of the line
            myNeighborTable->myNeighbors.push_back(neighbor);
        }
        
        for (int i=0; i<*numberOfNeighbors; i++){
            int neighbor;
            int numberOfSendNodes;
            meshFile >> buffer; //SEND
            meshFile >> buffer; //PE
            meshFile >> neighbor;
            meshFile >> numberOfSendNodes;
            std::getline(meshFile, buffer);//discard rest of the line
            //Assemble arrays of nodes to be sent
            for (int j=0; j<numberOfSendNodes; j++){
                int sendNode;
                meshFile >> sendNode;                
                myNeighborTable->myNeighbors[i].sendNodes.push_back(sendNode);
            }
            std::getline(meshFile, buffer);//discard rest of the line
        }
        
    }    

    void openfort14File(std::ifstream& meshFile, int domainID)
    {
        std::cerr << "Entering c subroutine openfort14File\n";
        std::string meshFileName = meshDir;
        meshFileName.append("/PE");
        std::stringstream buf;
        buf.width(4);
        buf.fill('0');
        buf << domainID;
        meshFileName.append(buf.str());
        meshFileName.append("/fort.14");
//        std::cerr << meshFileName << "\n";
        meshFile.open(meshFileName.c_str());
        if (!meshFile) {
            throw std::runtime_error("could not open fort.14 file "+meshFileName);
        }    
    }



    void readFort14Header(std::ifstream& meshFile, int *numberOfElements, int *numberOfPoints)
    {
        std::cerr << "Entering c subroutine readFort14Header\n";
        std::string buffer(1024, ' ');
        // discard first line, which only holds comments anyway
        std::getline(meshFile, buffer);

        meshFile >> *numberOfElements;
        meshFile >> *numberOfPoints;

        // discard remainder of line
        std::getline(meshFile, buffer);

//        std::cerr << "numberOfElements: " << *numberOfElements << "\n"
//                  << "numberOfPoints: " << *numberOfPoints << "\n";

        if (!meshFile.good()) {
            throw std::logic_error("could not read header");
        }
    }

    void readFort14Points(std::ifstream& meshFile, std::vector<FloatCoord<3> > *points, const int numberOfPoints)
    {
        std::cerr << "Entering c subroutine readFort14Points\n";

        std::string buffer(1024, ' ');

        for (int i = 0; i < numberOfPoints; ++i) {
            FloatCoord<3> p;
            int buf;

            meshFile >> buf;
            meshFile >> p[0];
            meshFile >> p[1];
            meshFile >> p[2];
            std::getline(meshFile, buffer);

            *points << p;
        }

        if (!meshFile.good()) {
            throw std::runtime_error("could not read points");
        }
    }


};




void runSimulation()
{
  
    Coord<2> dim(3, 2);
    std::size_t numCells = 100;
    double minDistance = 100;
    double quadrantSize = 400;
    quadrantDim = FloatCoord<2>(quadrantSize, quadrantSize);

    std::string prunedDirname("/home/zbyerly/adcirclgd/meshes/parallel_quarter_annular_v50_99");

    int steps = 100;

    SerialSimulator<ContainerCellType> sim(
        new ADCIRCInitializer(prunedDirname, steps));

    //TODO fix this stuff
    //    int ioPeriod = 1;
//    SiloWriter<ContainerCellType> *writer = new SiloWriter<ContainerCellType>("mesh", *ioPeriod);
//    writer->addSelectorForUnstructuredGrid(
//        &DomainCell::alive,
//        "DomainCell_wetdry");
//    sim.addWriter(writer);

    sim.run();
}

int main(int argc, char *argv[])
{
  //MPI_Init(&argc, &argv);
  //Typemaps::initializeMaps();

    runSimulation();

    //MPI_Finalize();
    return 0;
}

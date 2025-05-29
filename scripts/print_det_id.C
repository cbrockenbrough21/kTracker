#include <TSystem.h>
R__LOAD_LIBRARY(geom_svc)

int detName()
{
    // Initialize geometry service
    GeomSvc::UseDbSvc(true);
    GeomSvc* geom_svc = GeomSvc::instance();

    for (int detID = 0; detID <= 60; ++detID)
    {   
        std::string detName = geom_svc->getDetectorName(detID);
        cout << "Detector ID: " << detID << " => Name: " << detName  << "cell width: "<< geom_svc->getCellWidth (detID)  << endl;
    }   

    return 0;
}
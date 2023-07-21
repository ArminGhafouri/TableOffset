#include<Offset.hxx>
#include<TopExp_Explorer.hxx>
#include<Tools.hxx>
#include <STEPControl_Writer.hxx>
#include <Bspline.hxx>

///gpt
#include <iostream>
#include <vector>
#include <map>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>

#include <iostream>
#include <TopoDS_Face.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopExp.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>


using namespace std;

int main()
{
	try
	{

		Offset ofs(2);
	
		//ofs.PlotSurfaceOffsetOneByInterpolateCurve();

		ofs.PlotSurfaceOffsetOneByOccBsplineCurve();

		//ofs.PlotSurfaceOffsetOneByMyBspline();

		//ofs.PlotCurve(ofs.CreatCurveOffsetOneByInterpolate(), "Curve by interpolate", 10);

		//ofs.PlotSurfaceOffsetOneBySweepApp();

	}


	catch (const Standard_Failure& ex)
	{

		std::cout << ex.GetMessageString() << std::endl;


	}



	return 0;

}

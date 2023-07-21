#include"Offset.hxx"
#include<BRep_Tool.hxx>
#include<Tools.hxx>
#include<Bspline.hxx>




#include <vector>
#include<iostream>
#include <fstream>  
#include <string> 


#include <GeomFill_Coons.hxx> 
#include <GeomFill_Curved.hxx > 
#include <GeomFill_Stretch.hxx> 
#include <GeomConvert_CompBezierSurfacesToBSplineSurface.hxx> 
#include<Geom_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include<Geom2dAPI_InterCurveCurve.hxx>
#include <GeomAPI.hxx>
#include <gp_Pln.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomConvert.hxx>
#include <Geom_Plane.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <gp_Pnt2d.hxx>
#include <BRepBuilderAPI_MakeShape.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <NCollection_Array1.hxx>
#include <GeomFill_Sweep.hxx>
#include <GeomFill_SweepSectionGenerator.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <GeomFill_AppSweep.hxx>



#include<BrepOffsetAPI_MakeOffset.hxx>
#include <STEPControl_Writer.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Geometry.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepAlgoAPI_Fuse.hxx>


#include <GeomFill_BezierCurves.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepFill_CurveConstraint.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <Geom_BoundedSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomFill_Pipe.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom2dAPI_PointsToBSpline.hxx>
#include <TColGeom_Array1OfCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <Geom_Line.hxx>
#include <gp_Lin.hxx>
#include <GeomFill.hxx>
#include <TopoDS.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
#include <BRepPrimAPI_MakeSweep.hxx>
#include <GeomFill_SectionGenerator.hxx>
#include <GeomFill_AppSurf.hxx>
#include <GeomFill_Line.hxx>





Offset::Offset(int deg)
{
	theDegree_ = deg;
	/*theNumberOfSection_ = sec;
	theNumberOfPoint_ = num;*/
	//CalcTheU();

}

std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetOneByInterpolate() const
{

	std::vector<Handle(Geom_Curve)> curvetotal0;
	Handle(TColgp_HArray1OfPnt) Pts0 =
		new TColgp_HArray1OfPnt;

	std::vector <std::vector < gp_Pnt >> poinOff = ReadOffsetOne();



	for (int i = 0; i < poinOff.size(); i++)
	{

		Pts0->Resize(1, poinOff[i].size(), false);


		for (int j = 0; j < poinOff[i].size(); j++)
		{

			double x = poinOff[i][j].X();
			double y = poinOff[i][j].Y();
			double z = poinOff[i][j].Z();

			gp_Pnt pointfinall(x, y, z);

			Pts0->SetValue(j + 1, pointfinall /*poinOff[i][j].XYZ()*/);



			//std::cout <<"null ;"<< cur.IsNull();

		}

		/*GeomAPI_PointsToBSpline interpolator(Pts0);*/


		GeomAPI_Interpolate inter0(Pts0, false, 1e-6);
	
		inter0.Perform();
		Handle(Geom_Curve) cur0 = inter0.Curve();
		curvetotal0.push_back(cur0);


	}

	//std::cout << "*****null*****///////end0000;"  /*<< cur.IsNull() */ << std::endl;
	return curvetotal0;

}






void Offset::PlotSurfaceOffsetOneByGeomFill_BSplineCurves()const
{

	std::vector<TColgp_Array1OfPnt> points;
	std::vector<GeomAPI_PointsToBSpline> BSpoints;
	std::vector<Handle(Geom_BSplineCurve)> BScurve;
	//std::vector <std::vector < gp_Pnt >> poinOff = ReadOffsetOne();

	for (int i = 0; i < ReadOffsetOne().size(); i++)
	{


		TColgp_Array1OfPnt po /*= TColgp_Array1OfPnt(1 , vectorOfvector()[n].size() )*/;
		po.Resize(1, ReadOffsetOne()[i].size(), true);



		for (int j = 0; j < ReadOffsetOne()[i].size(); j++)
		{
			/*std::cout << "siZe vecofvec[I] =  " << vectorOfvector()[num].size() << std::endl;
			std::cout << "siZe vecofvec =  " << vectorOfvector().size() << std::endl;*/

			po.SetValue(j + 1, ReadOffsetOne()[i][j].XYZ());


		}

		/*std::cout << " // ********************* //  " << std::endl;
		std::cout << " //   *****************   //  " << std::endl;*/

		points.push_back(po);

	}


	/* std::cout << "  size points ======================>>>>>  " << points.size() << std::endl;
	 std::cout << "  size BScurve ======================>>>>>  " << BScurve.size() << std::endl;*/


	for (int r = 0; r < points.size(); r++)
	{

		/*    GeomAbs_C0,
			GeomAbs_G1,
			GeomAbs_C1,
			GeomAbs_G2,
			GeomAbs_C2,
			GeomAbs_C3,
			GeomAbs_CN*/

		GeomAPI_PointsToBSpline Bs;
		int DegMax = theDegree_;
		Bs = GeomAPI_PointsToBSpline(points[r], 0, DegMax /*5*/, GeomAbs_C1, 0.001);

		//BSpoints.push_back(Bs);
		Handle(Geom_BSplineCurve) cur;

		cur = Bs.Curve();
		//cur->IncreaseDegree(10);


		BScurve.push_back(cur);

	}


	//std::cout << "  size BScurve ===> " << BScurve.size() << std::endl;

	/*std::cout << points[3].Value(3).XYZ().X() << "  " << points[3].Value(3).XYZ().Y() << "  "
		<< points[3].Value(3).XYZ().Z() << std::endl;*/


	GeomFill_BSplineCurves fill;


	std::vector<TopoDS_Face> Faces;

	for (int i = 0; i < BScurve.size() - 1; i++)
	{

		//GeomFill_StretchStyle
		//GeomFill_CurvedStyle  
		//GeomFill_CoonsStyle


		fill.Init(BScurve[i], BScurve[i + 1], GeomFill_CoonsStyle);

		Handle(Geom_BSplineSurface) BSS = fill.Surface();
		//BSS->IncreaseDegree(6, 6);

		BRepBuilderAPI_MakeFace face(BSS, 5);
		Faces.push_back(face);


		/* else
		 {
			 std::cout << "ERROR" << std::endl;
		 }*/
	}


	TopoDS_Shape shape = Tools::JoinOCShape(Faces);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);

	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);
	BRepBuilderAPI_Transform build(shape, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();


	TopoDS_Shape FullShape = BRepAlgoAPI_Fuse(shape/*shape*/, halfshape).Shape();


	//Tools::PlotShapeTwo(FullShape, "PlotSurfaceOffsetOneByGeomFill_BSplineCurves");
	Tools::PlotShape(FullShape, "PlotSurfaceOffsetOneByGeomFill_BSplineCurves", 120, 120);


}



std::vector<std::vector<gp_Pnt>> Offset::ReadOffsetTwo() const
{
	//int result = -1;

	std::string theFileName_ = "geomet.out";
	std::fstream fin(theFileName_, std::ios::in);

	std::vector<double>XSectionTot;
	std::vector<double>NumberOfPoints;
	std::vector< std::vector<double>>Y_Coretotal;
	std::vector< std::vector<double>>Z_Coretotal;

	if (!fin.is_open())
	{
		std::exception e("Cannot Open Input .Out File");
		throw e;
	}

	std::string a;


	int NSection;
	char symm;
	double referenceDraft;

	fin >> NSection >> symm >> referenceDraft;
	bool symmetry;
	if (symm == 'T')
		symmetry = true;
	else if (symm == 'F')
		symmetry = false;

	//std::cout << NSection << "  " << symm << "  " << referenceDraft << std::endl;//
	// 
	// 
	// 
	//for (int i = 0; i < 3; i++)//
	//{
	//	getline(fin, a);
	//	//std::cout << a << std::endl;
	//}

	for (int i = 0; i < NSection; i++)
	{
		int nPoints, nGaps;
		double XSection;



		fin >> XSection >> nPoints >> nGaps;

		//thenGap_.push_back(nGaps);
		// 
		// 
		//std::cout << XSection << "  " << nPoints << "  " << nGaps << std::endl;/*'''''''''''''''''''''''''*/



		XSectionTot.push_back(XSection);
		NumberOfPoints.push_back(nPoints);


		// theGapindex
		std::vector<int> indeces = { 0 };
		//if (nGaps > 0) clear.indeces; //??
		for (int j = 0; j < nGaps; j++)
		{
			int index;
			fin >> index;
			indeces.push_back(index);
			//std::cout << index << "  ";
		}
		//std::cout << std::endl;
		//theGap_ .push_back(indeces);
		getline(fin, a);




		std::vector<double> yVec;
		std::vector<double> zVec;


		for (int j = 0; j < nPoints; j++)
		{
			double yValue;
			fin >> yValue;
			yVec.push_back(yValue);
			//std::cout << yValue << "  ";
		}

		Y_Coretotal.push_back(yVec); ///my
		//std::cout << std::endl;




		for (int j = 0; j < nPoints; j++)
		{
			double zValue;
			fin >> zValue;
			zVec.push_back(zValue);
			//std::cout << zValue << "  ";
		}
		//std::cout << std::endl;
		Z_Coretotal.push_back(zVec);


		//for (int j = 0; j < nPoints; j++)
		//{
		//	StripTheory_Point* Pt = new StripTheory_Point
		//	(
		//		XSection,
		//		yVec[j],
		//		zVec[j]
		//	);
		//	thePoints_.push_back(Pt);
		// 
		//}




	}

	//for (int i = 0; i < XSectionTot.size(); i++)
	//{

	//	std::cout << " SECTION{" << i << "} ==> " << XSectionTot[i] << std::endl;

	//	for (int j = 0; j < Y_Coretotal[i].size(); j++)
	//	{

	//		std::cout << /*XSectionTot[i] << "  " <<*/ Y_Coretotal[i][j] << "   ,   " << Z_Coretotal[i][j]  /*<< Z_Coretotal[i][j] <<*/ << std::endl;

	//	}
	//	std::cout << std::endl;

	//	//std::cout << "//////////////////// " << std::endl;

	//}




	std::vector< std::vector<gp_Pnt> > pointss0;


	for (int i = 0; i < NSection; i++)
	{


		std::vector<gp_Pnt> points0;
		double xdummy = XSectionTot[i];


		for (int j = 0; j < Y_Coretotal[i].size(); j++)
		{

			double Ysec = Y_Coretotal[i][j];
			double zdummy = Z_Coretotal[i][j];


			gp_Pnt point1(xdummy, Ysec, zdummy);

			points0.push_back(point1);


		}


		pointss0.push_back(points0);

	}

	return pointss0;
}


///**********************************************************************************////
//Creat curves from tableOffset 3

std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetTwoByInterpolate() const
{
	/*std::cout << "//////////////////// 2  " << std::endl;*/
	//std::vector<Handle(TColgp_HArray1OfPnt)> PtsT;
	std::vector<std::vector<gp_Pnt>> poinOff = ReadOffsetTwo();
	std::vector<Handle(Geom_Curve)> curvetotal0;

	Handle(TColgp_HArray1OfPnt) Pts0 = new TColgp_HArray1OfPnt;

	//Pts0->Resize(1, vectorOfvector2()[0].size()/*vectorOfvector2()[i].size()*/, false);


	for (int i = 0; i < poinOff.size(); i++)
	{


		Pts0->Resize(1, poinOff[i].size(), false);
		//std::cout << "////////////////////66 " << std::endl;



		for (int j = 0; j < poinOff[i].size()/*vectorOfvector2()[i].size()*/; j++)
		{
			Pts0->SetValue(j + 1, poinOff[i][j]/*.XYZ()*/);

		}

		GeomAPI_Interpolate inter0(Pts0, false, 0.01);
		inter0.Perform();
		Handle(Geom_Curve) cur0 = inter0.Curve();
		curvetotal0.push_back(cur0);
	}

	//std::cout << "*****null*****///////end0000;"  /*<< cur.IsNull() */ << std::endl;
	return curvetotal0;
}


///***************************************************************************/////
//Plot Curves from CreatCurve2() or tableoffset 3;



void Offset::PlotSurfaceOffsetTwoByGeomFill_BSplineCurves() const
{

	std::vector<TColgp_Array1OfPnt> points;
	std::vector<GeomAPI_PointsToBSpline> BSpoints;
	std::vector<Handle(Geom_BSplineCurve)> BScurve;


	for (int i = 0; i < ReadOffsetTwo().size(); i++)
	{

		//int num = i;
		TColgp_Array1OfPnt po;
		po.Resize(1, ReadOffsetTwo()[i].size(), true);



		for (int j = 0; j < ReadOffsetTwo()[i].size(); j++)
		{
			/*std::cout << "siZe vecofvec[I] =  " << vectorOfvector()[num].size() << std::endl;
			std::cout << "siZe vecofvec =  " << vectorOfvector().size() << std::endl;*/

			po.SetValue(j + 1, ReadOffsetTwo()[i][j].XYZ());

			// std::cout << vectorOfvector()[i][j].XYZ().X() << "  " << vectorOfvector()[i][j].XYZ().Y() << "  "
			//<< vectorOfvector()[i][j].XYZ().Z() << std::endl;

		   /* std::cout << points[i].Value(j + 1).XYZ().X() << "  " << points[i].Value(j + 1).XYZ().Y() << "  "
				<< points[i].Value(j+1).XYZ().Z() << std::endl;*/
				//points.push_back(po);

				//std::cout << points[i].Value(j+1).XYZ().X() << "  " << points[i].Value(j+1).XYZ().Y() << "  "
				   // << points[i].Value(j+1).XYZ().Z() << std::endl;

			std::cout << ReadOffsetTwo()[i][j].XYZ().X() << "  " << ReadOffsetTwo()[i][j].XYZ().Y() << "  "
				<< ReadOffsetTwo()[i][j].XYZ().Z() << std::endl;

		}

		std::cout << " // ********************* //  " << std::endl;
		std::cout << " //   *****************   //  " << std::endl;

		points.push_back(po);




	}


	/* std::cout << "  size points ======================>>>>>  " << points.size() << std::endl;
	 std::cout << "  size BScurve ======================>>>>>  " << BScurve.size() << std::endl;*/


	for (int r = 0; r < points.size(); r++)
	{

		/*    GeomAbs_C0,
			GeomAbs_G1,
			GeomAbs_C1,
			GeomAbs_G2,
			GeomAbs_C2,
			GeomAbs_C3,
			GeomAbs_CN*/

		GeomAPI_PointsToBSpline Bs;
		Bs = GeomAPI_PointsToBSpline(points[r], 0, 2, GeomAbs_C1, 0.001);

		//BSpoints.push_back(Bs);
		Handle(Geom_BSplineCurve) cur;
		cur = Bs.Curve();


		BScurve.push_back(cur);

	}


	//std::cout << "  size BScurve====> " << BScurve.size() << std::endl;

	std::cout << points[3].Value(3).XYZ().X() << "  " << points[3].Value(3).XYZ().Y() << "  "
		<< points[3].Value(3).XYZ().Z() << std::endl;



	GeomFill_BSplineCurves fill;
	//GeomFill_BSplineCurves fill2;
	std::vector<TopoDS_Face> Faces;
	//std::vector<Handle(Geom_BSplineSurface)> surfaces;

	for (int i = 0; i < BScurve.size() - 1; i++)
	{

		//GeomFill_StretchStyle
		//GeomFill_CurvedStyle  
		//GeomFill_CoonsStyle


		fill.Init(BScurve[i], BScurve[i + 1], GeomFill_CoonsStyle);

		Handle(Geom_BSplineSurface) BSS = fill.Surface();
		//BSS->IncreaseDegree(3, 5);
		BRepBuilderAPI_MakeFace face(BSS, 5);
		Faces.push_back(face);


		/* else
		 {
			 std::cout << "ERROR" << std::endl;
		 }*/

	}


	TopoDS_Shape shape = Tools::JoinOCShape(Faces);



	gp_Pnt po0(ReadOffsetTwo()[0][0].X(), ReadOffsetTwo()[0][0].Y(), ReadOffsetTwo()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);

	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);
	BRepBuilderAPI_Transform build(shape, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();


	TopoDS_Shape FullShape = BRepAlgoAPI_Fuse(shape/*shape*/, halfshape).Shape();
	//Tools::PlotShapeTwo(shape, "PlotSurfaceOffsetTwoByGeomFill_BSplineCurves");
	Tools::PlotShape(FullShape, "PlotSurfaceOffsetTwoByGeomFill_BSplineCurves", 100, 100);


}

//**********************************************************************************///
//Creat Bspline curves from tableOffset 3

std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetTwoBySpline() const
{
	std::vector<std::vector<gp_Pnt>> poinOff = ReadOffsetTwo();
	std::vector<Handle(Geom_Curve)> curvetotal;


	for (int i = 0; i < poinOff.size(); i++)
	{
		std::vector<gp_Pnt> Poin;

		for (int j = 0; j < poinOff[i].size(); j++)
		{

			gp_Pnt pnt = poinOff[i][j];
			Poin.push_back(pnt);

		}


		Bspline bs(theDegree_, Poin);
		Handle(Geom_Curve) cur = bs.CreateOCCurve();

		curvetotal.push_back(cur);

	}

	return  curvetotal;

}

//************************************************************************************////
///creat surface by correction tableOffset 2 
//Removing extra zeros in the sections , and set y = 0.0 at the beginning of sections

void Offset::PlotSurfaceOffsetOneByInterpolateCurve() const
{

	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveOffsetOneByInterpolate/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{

		col.Append(CreatCurveOffsetOneByInterpolate/*CreatCurveOffsetOneBySpline*/()[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(2, 2, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree() , anAlgo.VDegree());


	//BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	BRepBuilderAPI_MakeFace face(aRes,1 ,2  ,1 ,2 , 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);

	//TopoDS_Shell shell = TopoDS::Shell(ship);
	/*TopoDS_Shape shapeofship = ship.Complemented();
	TopoDS_Solid solid = TopoDS::Solid(shapeofship);*/
	//opoDS_Shape FullShape = BRepAlgoAPI_Fuse(airplane_compound/*shape*/, halfshape).Shape();


   //Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceOffsetOneByInterpolateCurve", 100, 100);

}

//************************************************************************************///


std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetOneBySpline() const
{

	std::vector <std::vector < gp_Pnt >> poinOff = ReadOffsetOne();

	std::vector<Handle(Geom_Curve)> curvetotal;

	for (int i = 0; i < poinOff.size(); i++)
	{
		std::vector<gp_Pnt> Poin;

		for (int j = 0; j < poinOff[i].size(); j++)
		{

			gp_Pnt pnt = poinOff[i][j];
			Poin.push_back(pnt);

		}

		Bspline bs(theDegree_, Poin);
		Handle(Geom_Curve) cur = bs.CreateOCCurve();

		curvetotal.push_back(cur);

	}


	return  curvetotal;

}




void Offset::PlotSurfaceOffsetTwoByInterpolateCurve() const
{
	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveOffsetTwoByInterpolate().size(); i++)
	{

		col.Append(CreatCurveOffsetTwoByInterpolate()[i]);

	}


	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 0.0000001, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, 5, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 0.1);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetTwo()[0][0].X(), ReadOffsetTwo()[0][0].Y(), ReadOffsetTwo()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	//Tools::PlotShapeTwo(ship, "Shape BY curve" /*, 10 , 20*/);
	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceOffsetTwoByInterpolateCurve", 50, 50);

}



//*********************************************************************************//
// Get Correction points from tableOffset 2

std::vector<std::vector<gp_Pnt>> Offset::ReadOffsetOne() const
{
	std::fstream fin("Parent+49.txt", std::ios::in);

	if (!fin.is_open())
	{
		std::exception e("file not open");
		throw e;

	}


	std::string a;
	int numberofWaterline, numberofeStation;
	double  length, beam, draft;

	fin >> numberofWaterline;
	getline(fin, a);
	fin >> numberofeStation;

	/*std::cout << "numberofWaterline = "<< numberofWaterline << "     :    " <<"numberofeStation = "
	<< numberofeStation << std::endl;*/


	for (int i = 0; i < 6; i++)
	{

		getline(fin, a);

	}


	fin >> length;
	getline(fin, a);
	fin >> beam;
	getline(fin, a);
	fin >> draft;
	getline(fin, a);
	std::cout << length << " ; " << beam << " ; " << draft << std::endl;

	std::vector<double> Z_Core;
	std::vector<double>X_Core;
	std::vector< std::vector<double>>Y_Coretotal;





	for (int i = 0; i < numberofWaterline; i++)
	{
		double z;
		fin >> z;


		Z_Core.push_back(z);
		//std::cout << Z_Core[i] << std::endl;

	}
	// std::cout << Z_Core[15] << std::endl;



	for (int i = 0; i < 2; i++)
	{

		getline(fin, a);

	}
	//std::cout << a << std::endl;



	for (int i = 0; i < numberofeStation; i++)
	{

		double x;
		fin >> x;
		X_Core.push_back(x);

		//std::cout << x << std::endl;

		std::vector<double>Y_Core;

		for (int j = 0; j < numberofWaterline; j++)
		{
			double y;
			fin >> y;


			Y_Core.push_back(y);


		}

		Y_Coretotal.push_back(Y_Core);

		//std::cout << Y_Coretotal[i][1] << std::endl;
		getline(fin, a);

	}



	std::vector< std::vector<gp_Pnt> > pointss;


	for (int i = 0; i < numberofeStation; i++)
	{


		std::vector<gp_Pnt> points;
		double xdummy = X_Core[i];


		for (int j = 0; j < numberofWaterline; j++)
		{

			double Ysec = Y_Coretotal[i][j];
			double zdummy = Z_Core[j];

			/*if (j == 0)
			{
				Ysec = 0.0;
			}*/

			gp_Pnt point1(xdummy, Ysec, zdummy);

			points.push_back(point1);


		}


		pointss.push_back(points);

	}


	std::vector< std::vector<gp_Pnt> > pointss2;

	for (int i = 0; i < pointss.size(); i++)
	{


		std::vector<gp_Pnt> points0;
		double xdummy = X_Core[i];


		for (int j = 0; j < pointss[i].size() /*-1*/; j++)
		{

			double y1 = Y_Coretotal[i][j];


			if (std::abs(y1) != 0.0 && j != 0)
			{


				double Ysec000 = Y_Coretotal[i][j - 1];
				double zdummy = Z_Core[j - 1];

				gp_Pnt point00(xdummy, Ysec000, zdummy);
				points0.push_back(point00);

			}


		}

		int lastpoint = pointss[i].size() - 1;
		double Ysec00 = Y_Coretotal[i][lastpoint];
		double zdummy = Z_Core[lastpoint];
		gp_Pnt point00(xdummy, Ysec00, zdummy);
		points0.push_back(point00);
		pointss2.push_back(points0);
	}


	std::vector< std::vector<gp_Pnt> > pointssFinall;


	for (int i = 0; i < pointss2.size(); i++)
	{


		std::vector<gp_Pnt> points;
		double xdummy = X_Core[i];


		for (int j = 0; j < pointss2[i].size(); j++)
		{

			double Ysec0 = pointss2[i][j].Y();
			double zdummy0 = pointss2[i][j].Z();

			if (j == 0)
			{
				Ysec0 = 0.0;
			}

			gp_Pnt point1(xdummy, Ysec0, zdummy0);

			points.push_back(point1);


		}


		pointssFinall.push_back(points);

	}

	for (int i = 0; i < pointss2.size(); i++)
	{

		//std::cout << "*************** NEW SEC *************" << std::endl;
		//std::cout << "*********************" << std::endl;

		//for (int j = 0; j < pointssFinall[i].size(); j++)
		//{
		//	//std::cout << " size two (" << i << ") = " << pointss2[i].size() << std::endl;
		//	std::cout << " Section [" << i << "] = " << pointssFinall[i][j].X()
		//		<< "  " << pointssFinall[i][j].Y() << "  " << pointssFinall[i][j].Z() << std::endl;
		//}

	}

	return pointssFinall;
}




TopoDS_Shape Offset::CreatShapeFromTableOffset1()const
{
	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveOffsetOneByInterpolate().size(); i++)
	{

		col.Append(CreatCurveOffsetOneByInterpolate()[i]);

	}


	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 0.0001, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, aMaxDeg, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	return ship;
}




TopoDS_Shape Offset::CreatShapeFromTableOffset2()const
{

	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveOffsetTwoByInterpolate() /*or -> CreatCurve3() */.size(); i++)
	{

		col.Append(CreatCurveOffsetTwoByInterpolate()[i] /*or -> CreatCurve3()[i] */);

	}


	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 0.0001, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, 2, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	return ship;
}












//void Offset::PlotSurfaceOffsetOneByBsOCCCurve() const
//{
//	NCollection_List < Handle(Geom_Curve)> col;
//
//	for (int i = 0; i < CreatCurveOffsetOneBySpline/*CreatCurveOffsetOneBySpline*/().size(); i++)
//	{
//
//		col.Append(CreatCurveOffsetOneBySpline/*CreatCurveOffsetOneBySpline*/()[i]);
//
//	}
//
//	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
//	{*/
//	std::cout << "Ncollection_list size: " << col.Size();
//	//populate section generator
//
//	//GeomFill_SweepSectionGenerator aSecGenerator;
//	GeomFill_SectionGenerator aSecGenerator;
//
//	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
//	{
//		const Handle(Geom_Curve)& aCurve = anIt.Value();
//		aSecGenerator.AddCurve(aCurve);
//	}
//	aSecGenerator.Perform(Precision::PConfusion());
//	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());
//
//	//parameters
//	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
//	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);
//
//	//algorithm
//	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
//	GeomFill_AppSurf anAlgo(1, 5, aTol3d, aTol2d, aNbIt, false);
//
//	//anAlgo.Perform(aLine, aSecGenerator);
//	anAlgo.Perform(aLine, aSecGenerator);
//	Handle(Geom_Surface) aRes;
//
//	/*if (!anAlgo.IsDone())
//	{
//		return aRes;
//	}*/
//
//	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
//		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
//		anAlgo.UDegree(), anAlgo.VDegree());
//
//
//	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
//	TopoDS_Shape shape = face.Shape();
//
//
//	BRep_Builder builder;
//	TopoDS_Compound ship;
//	builder.MakeCompound(ship);
//	builder.Add(ship, shape);
//
//
//	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
//	gp_Dir dir0(0.0, 1.0, 0.0);
//	gp_Dir dir1(0.0, 0.0, 1.0);
//	gp_Ax2 ax0(po0, dir0, dir1);
//	gp_Trsf trsf2;
//	trsf2.SetMirror(ax0);
//
//	BRepBuilderAPI_Transform build(ship, trsf2, true);
//	TopoDS_Shape halfshape = build.Shape();
//
//	builder.Add(ship, halfshape);
//
//	//TopoDS_Shell shell = TopoDS::Shell(ship);
//	/*TopoDS_Shape shapeofship = ship.Complemented();
//	TopoDS_Solid solid = TopoDS::Solid(shapeofship);*/
//	//opoDS_Shape FullShape = BRepAlgoAPI_Fuse(airplane_compound/*shape*/, halfshape).Shape();
//
//
//   //Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
//	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceOffsetOneByBsplineCurveTest", 20 , 20 );
//
//}


void Offset::PlotSurfaceOffsetOneBySweepApp() const
{
	NCollection_List < Handle(Geom_Curve)> col;
	for (int i = 0; i < CreatCurveOffsetOneByInterpolate/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{
		col.Append(CreatCurveOffsetOneByInterpolate/*CreatCurveOffsetOneBySpline*/()[i]);

	}
	

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/

	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator


	GeomFill_SweepSectionGenerator aSecGenerator;


	for (int j = 0; j < CreatCurveOffsetOneByInterpolate().size() - 1; j++)
	{
		Handle(Geom_Curve) apath = CreatCurvePathOffsetOneByInterpolateBottomeVec()[0];
		Handle(Geom_Curve) aCurvfirst0 = CreatCurveOffsetOneByInterpolate()[j];
		Handle(Geom_Curve) aCurvelast = CreatCurveOffsetOneByInterpolate()[j + 1]; 


	/*	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
		{
			Handle(Geom_Curve) aCurve = anIt.Value();*/


			/*Handle(Geom_Curve) aCurvfirst = CreatCurveOffsetOneByInterpolate()[j];
			Handle(Geom_Curve) aCurvelast = CreatCurveOffsetOneByInterpolate()[j+1];*/
				

	
			aSecGenerator.Init(apath/*, aCurve*/, aCurvfirst0, aCurvelast);

		
		}
	/*}*/

	aSecGenerator.Perform(Precision::PIntersection()/*::PConfusion()*/);
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0;
	Standard_Real aTol3d = 1e-6, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);


	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes/*, 1, 3, 1, 3 */, 0.00001);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	Tools::PlotShape(ship/*shell.Complemented()*/, "SWEEpPlotSurfaceOffsetOneByBsplineCurveTest", 20, 20);
	//Tools::PlotShapeTwo(ship/*shell.Complemented()*/, "SWEEpPlotSurfaceOffsetOneByBsplineCurveTest");


}













std::vector<Handle(Geom_Curve)> Offset::CreatCurveTheBottomOftheShip() const
{
	std::vector<std::vector<gp_Pnt>> points = ReadOffsetOne();
	//std::vector<Handle(Geom_Curve)> curvetotal;

	/*for (int i = 0; i < points.size(); i++)
	{


		for (int j = 0; j < points[i].size(); j++)
		{
			if (j == 0 || j == 1)
			{
			

			}

		}



	}*/
	
	Handle(TColgp_HArray1OfPnt) Pts =
		new TColgp_HArray1OfPnt;


	std::vector<Handle(Geom_Curve)> curvetotal;
	

	for (int i = 2; i < points.size()-3 ; i++)
	{
		Pts->Resize(1, 2, false);

		for (int j = 0; j < points[i].size() ; j++)
		{

			if (j == 0 || j == 1)
			{
				double X = points[i][j].X();
				double Y = points[i][j].Y();
				double Z = points[i][j].Z();;
				gp_Pnt points2(X, Y, Z);
				Pts->SetValue(j + 1, points2);

				//Pts->SetValue(j + 1, points[i][j]);
			}
			

		}

		GeomAPI_Interpolate inter(Pts, false, 0.0001);
		inter.Perform();
		Handle(Geom_Curve) cur = inter.Curve();

		curvetotal.push_back(cur);

	}


	return  curvetotal;

}



void Offset::PlotCurveTheBottomOftheShip(std::string filename, int n)
{
	std::fstream My_File(filename, std::ios::out);
	if (!My_File.is_open())
	{
		std::cout << "file not open";
		return;
	}


	My_File << "VARIABLES = X Y Z" << std::endl;
	// My_File << "ZONE T = Curve" << std::endl;



	auto curve0 = CreatCurveTheBottomOftheShip();


	for (int j = 0; j < curve0.size(); j++)
	{

		My_File << "ZONE T = Curve" << std::endl;
		double du = (curve0[j]->LastParameter() - curve0[j]->FirstParameter()) / n;
		/* std::cout << curve0[j]->FirstParameter() << std::endl;
		 std::cout << curve0[j]->LastParameter();*/

		for (int i = 0; i <= n; i++)
		{
			double u = (i * du) + curve0[j]->FirstParameter();
			gp_Pnt p = curve0[j]->Value(u);


			My_File << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;

		}




	}

	My_File.close();

}



void Offset::PlotSurfaceTheBottomOftheShipByInterpolateCurve() const
{

	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveTheBottomOftheShip/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{

		col.Append(CreatCurveTheBottomOftheShip/*CreatCurveOffsetOneBySpline*/()[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(0, 5, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);

	//TopoDS_Shell shell = TopoDS::Shell(ship);
	/*TopoDS_Shape shapeofship = ship.Complemented();
	TopoDS_Solid solid = TopoDS::Solid(shapeofship);*/
	//opoDS_Shape FullShape = BRepAlgoAPI_Fuse(airplane_compound/*shape*/, halfshape).Shape();


   //Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceONEBottomeShip", 100, 100);



}




std::vector<Handle(Geom_Curve)> Offset::CreatCurveTheCrustOftheShip() const
{
	std::vector<std::vector<gp_Pnt>> points = ReadOffsetOne();

	Handle(TColgp_HArray1OfPnt) Pts =
		new TColgp_HArray1OfPnt;


	std::vector<Handle(Geom_Curve)> curvetotal;


	for (int i = 0; i < points.size(); i++)
	{
		
		Pts->Resize(1, points[i].size() - 1 , false);

		for (int j = 1 ; j < points[i].size(); j++)
		{

			Pts->SetValue( j  , points[i][j]);
			/*std::cout << "null ;" << cur.IsNull(); */
		}

		GeomAPI_Interpolate inter(Pts, false, 0.01);
		inter.Perform();
		Handle(Geom_Curve) cur = inter.Curve();


		curvetotal.push_back(cur);
	}

	return curvetotal;
 }





void Offset::PlotSurfaceTheCrustOftheShipByInterpolateCurve() const
{
	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveTheCrustOftheShip/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{

		col.Append(CreatCurveTheCrustOftheShip/*CreatCurveOffsetOneBySpline*/()[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(0, 5, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);

	//TopoDS_Shell shell = TopoDS::Shell(ship);
	/*TopoDS_Shape shapeofship = ship.Complemented();
	TopoDS_Solid solid = TopoDS::Solid(shapeofship);*/
	//opoDS_Shape FullShape = BRepAlgoAPI_Fuse(airplane_compound/*shape*/, halfshape).Shape();


   //Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceCRUST", 100, 100);



}

void Offset::PlotCurveTheCrustOftheShip(std::string filename, int n)
{

	std::fstream My_File(filename, std::ios::out);
	if (!My_File.is_open())
	{
		std::cout << "file not open";
		return;
	}


	My_File << "VARIABLES = X Y Z" << std::endl;
	// My_File << "ZONE T = Curve" << std::endl;



	auto curve0 = CreatCurveTheCrustOftheShip();


	for (int j = 0; j < curve0.size(); j++)
	{

		My_File << "ZONE T = Curve" << std::endl;
		double du = (curve0[j]->LastParameter() - curve0[j]->FirstParameter()) / n;
		/* std::cout << curve0[j]->FirstParameter() << std::endl;
		 std::cout << curve0[j]->LastParameter();*/

		for (int i = 0; i <= n; i++)
		{
			double u = (i * du) + curve0[j]->FirstParameter();
			gp_Pnt p = curve0[j]->Value(u);


			My_File << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;

		}




	}

	My_File.close();




}


TopoDS_Shape Offset::CreatShapeFromTableOffset1Bottome() const
{
	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveTheBottomOftheShip/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{

		col.Append(CreatCurveTheBottomOftheShip/*CreatCurveOffsetOneBySpline*/()[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(0, 5, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);



    return ship;
}








TopoDS_Shape Offset::CreatShapeFromTableOffset1Crust() const
{
	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveTheCrustOftheShip/*CreatCurveOffsetOneBySpline*/().size(); i++)
	{

		col.Append(CreatCurveTheCrustOftheShip/*CreatCurveOffsetOneBySpline*/()[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size();
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(0, 5, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Perform(aLine, aSecGenerator);
	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);
	return ship;
}

void Offset::ShapeOfShip() const
{

	BRep_Builder builder;
	TopoDS_Compound ship_compound;
	builder.MakeCompound(ship_compound);
	builder.Add(ship_compound, CreatShapeFromTableOffset1Bottome());
	builder.Add(ship_compound, CreatShapeFromTableOffset1Crust());
	Tools::PlotShape(ship_compound, "Creatship",10,10);
}


void Offset::PlotSurfaceOffsetOneTESTBottome() const
{
	std::vector<TColgp_Array1OfPnt> points;
	std::vector<GeomAPI_PointsToBSpline> BSpoints;
	std::vector<Handle(Geom_BSplineCurve)> BScurve;
	std::vector<std::vector<gp_Pnt>> pointof = ReadOffsetOne();

	for (int i = 0; i < pointof.size(); i++)
	{

		TColgp_Array1OfPnt po /*= TColgp_Array1OfPnt(1 , vectorOfvector()[n].size() )*/;
		po.Resize(1, 2 , true);

		for (int j = 0; j < pointof[i].size(); j++)
		{
			if (j == 0 || j == 1)
			{
				

				/*double X = pointof[i][j].X();
				double Y = pointof[i][j].Y();
				double Z = 0.0;
				gp_Pnt points2(X, Y, Z);
				
				po.SetValue(j + 1, points2);*/

				po.SetValue(j + 1, pointof[i][j]/*.XYZ()*/);
			}

		}

		std::cout << " // ********************* //  " << std::endl;
		std::cout << " //   *****************   //  " << std::endl;

		points.push_back(po);

	}



	for (int num = 0; num < points.size(); num++)
	{

		/*    GeomAbs_C0,
			GeomAbs_G1,
			GeomAbs_C1,
			GeomAbs_G2,
			GeomAbs_C2,
			GeomAbs_C3,
			GeomAbs_CN*/

		GeomAPI_PointsToBSpline Bs;
		int DegMax = theDegree_;
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_C0, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_G1, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_G2, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_C1, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_C2, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 /*5*/, GeomAbs_C3, 0.00001);
		Bs = GeomAPI_PointsToBSpline(points[num], 1, 3 , GeomAbs_C3, 0.00001);

		//BSpoints.push_back(Bs);
		Handle(Geom_BSplineCurve) cur;

		cur = Bs.Curve();
		//cur->IncreaseDegree(10);
		BScurve.push_back(cur);

	}


	NCollection_List < Handle(Geom_BSplineCurve)> col;
	//NCollection_TListIterator < Handle(Geom_BSplineCurve)> col0;

	for (int i = 0; i < BScurve.size(); i++)
	{

		col.Append(BScurve[i]);

	}

	
	std::cout << "Ncollection_list size: " << col.Size() << std::endl;
	//populate section generator

	//GeomFill_SweepSectionGenerator aSecGenerator;
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision/*::Angular()*/::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
	Standard_Real aTol3d = 1e-6, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, 3, aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Init(1 , 3, aTol3d, aTol2d, aNbIt, false);
	//anAlgo.Perform(aLine, aSecGenerator);

	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
	anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
	anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1e-6);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	//Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
	Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceOffsetOneByInterTestBottome", 10, 10);


}




std::vector<Handle(Geom_Curve)> Offset::CreatCurvePathOffsetOneByInterpolateBottomeVec() const
{

	std::vector<Handle(Geom_Curve)> curvetotal0;
	Handle(TColgp_HArray1OfPnt) Pts0 =
		new TColgp_HArray1OfPnt;

	std::vector <std::vector < gp_Pnt >> poinOff = ReadOffsetOne();

	for (int i = 0; i < poinOff.size(); i++)
	{

		Pts0->Resize(1, 2/*poinOff[i].size()*/, false);


		//for (int j = 0; j < poinOff[i].size(); j++)
		//{
		int num = 1;
			Pts0->SetValue(/*j +*/ num, poinOff[i][0]/*.XYZ()*/);

		
		//}



		num++;

	}
	GeomAPI_Interpolate inter0(Pts0, false, 1e-6);
	//inter0.ClearTangents();
	inter0.Perform();
	Handle(Geom_Curve) cur0 = inter0.Curve();
	curvetotal0.push_back(cur0);


	//std::cout << "*****null*****///////end0000;" << curvetotal0.size() << std::endl;
	return curvetotal0;



}

void Offset::PlotCurve(std::vector<Handle(Geom_Curve)> curve , std::string filename , int n)
{

	std::fstream My_File(filename + ".plt", std::ios::out);
	if (!My_File.is_open())
	{
		std::cout << "file not open";
		return;
	}


	My_File << "VARIABLES = X Y Z" << std::endl;
	// My_File << "ZONE T = Curve" << std::endl;

	auto curve0 = curve;


	for (int j = 0; j < curve0.size(); j++)
	{

		My_File << "ZONE T = Curve" << std::endl;
		double du = (curve0[j]->LastParameter() - curve0[j]->FirstParameter()) / n;
		/* std::cout << curve0[j]->FirstParameter() << std::endl;
		 std::cout << curve0[j]->LastParameter();*/

		for (int i = 0; i <= n; i++)
		{
			double u = (i * du) + curve0[j]->FirstParameter();
			gp_Pnt p = curve0[j]->Value(u);


			My_File << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;

		}


	}

	My_File.close();


}





void Offset::PlotSurfaceOffsetOneByOccBsplineCurve() const
{

	std::vector<TColgp_Array1OfPnt> points;
	std::vector<GeomAPI_PointsToBSpline> BSpoints;
	std::vector<Handle(Geom_BSplineCurve)> BScurve;


	for (int i = 0; i < ReadOffsetOne().size(); i++)
	{
	
		TColgp_Array1OfPnt po /*= TColgp_Array1OfPnt(1 , vectorOfvector()[n].size() )*/;
		po.Resize(1, ReadOffsetOne()[i].size(), true);



		for (int j = 0; j < ReadOffsetOne()[i].size(); j++)
		{
			/*std::cout << "siZe vecofvec[I] =  " << vectorOfvector()[num].size() << std::endl;
			std::cout << "siZe vecofvec =  " << vectorOfvector().size() << std::endl;*/

			po.SetValue(j + 1, ReadOffsetOne()[i][j].XYZ());


		}

		std::cout << " // ********************* //  " << std::endl;
		std::cout << " //   *****************   //  " << std::endl;

		points.push_back(po);

	}



	for (int num = 0; num < points.size(); num++)
	{

		  /*GeomAbs_C0,
			GeomAbs_G1,
			GeomAbs_C1,
			GeomAbs_G2,
			GeomAbs_C2,
			GeomAbs_C3,
			GeomAbs_CN*/

		GeomAPI_PointsToBSpline Bs; 
	
		//int DegMax = theDegree_;
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C0, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5, GeomAbs_G1, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_G2, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C1, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C2, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C3, 0.00001);

		Bs = GeomAPI_PointsToBSpline(points[num], theDegree_ /*, 5 */, GeomAbs_C0, 1e-6);
		
		const Handle(Geom_BSplineCurve)& cur = Bs.Curve();
		
		BScurve.push_back(cur);

	}


	NCollection_List < Handle(Geom_BSplineCurve)> col;
	//NCollection_TListIterator < Handle(Geom_BSplineCurve)> col0;

	for (int i = 0; i < BScurve.size(); i++)
	{

		col.Append(BScurve[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size() << std::endl;
	//populate section generator

	
	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::Approximation());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, 2  , aTol3d, aTol2d, aNbIt, false);

	//anAlgo.Init(1 , 3, aTol3d, aTol2d, aNbIt, false);
	//anAlgo.Perform(aLine, aSecGenerator);

	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 0.5);
	TopoDS_Shape shape = face.Shape();


	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, true);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	//Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
	Tools::PlotShape(ship, "PlotSurfaceOffsetOneByOccBsplineCurve2", 20, 20);

}


void Offset::PlotSurfaceOffsetOneByMyBspline() const
{


NCollection_List < Handle(Geom_Curve)> col;

for (int i = 0; i < CreatCurveOffsetOneBySpline/*CreatCurveOffsetOneBySpline*/().size(); i++)
{

	col.Append(CreatCurveOffsetOneBySpline/*CreatCurveOffsetOneBySpline*/()[i]);

}

/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
{*/
std::cout << "Ncollection_list size: " << col.Size();
//populate section generator

//GeomFill_SweepSectionGenerator aSecGenerator;
GeomFill_SectionGenerator aSecGenerator;

for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(/*theSections*/col); anIt.More(); anIt.Next())
{
	const Handle(Geom_Curve)& aCurve = anIt.Value();
	aSecGenerator.AddCurve(aCurve);
}
aSecGenerator.Perform(Precision::PConfusion());
Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

//parameters
const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0 /*0*/;
Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);

//algorithm
//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
GeomFill_AppSurf anAlgo(1, 5, aTol3d, aTol2d, aNbIt, false);

//anAlgo.Perform(aLine, aSecGenerator);
anAlgo.Perform(aLine, aSecGenerator);
Handle(Geom_Surface) aRes;

/*if (!anAlgo.IsDone())
{
	return aRes;
}*/

aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
	anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
	anAlgo.UDegree(), anAlgo.VDegree());


BRepBuilderAPI_MakeFace face(aRes, 1e-6);
TopoDS_Shape shape = face.Shape();


BRep_Builder builder;
TopoDS_Compound ship;
builder.MakeCompound(ship);
builder.Add(ship, shape);


gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
gp_Dir dir0(0.0, 1.0, 0.0);
gp_Dir dir1(0.0, 0.0, 1.0);
gp_Ax2 ax0(po0, dir0, dir1);
gp_Trsf trsf2;
trsf2.SetMirror(ax0);

BRepBuilderAPI_Transform build(ship, trsf2, true);
TopoDS_Shape halfshape = build.Shape();

builder.Add(ship, halfshape);

//TopoDS_Shell shell = TopoDS::Shell(ship);
/*TopoDS_Shape shapeofship = ship.Complemented();
TopoDS_Solid solid = TopoDS::Solid(shapeofship);*/
//opoDS_Shape FullShape = BRepAlgoAPI_Fuse(airplane_compound/*shape*/, halfshape).Shape();


//Tools::PlotShapeTwo(ship, "PlotSurfaceOffsetOneByInterpolateCurve" /*, 10 , 20*/);
Tools::PlotShape(ship/*shell.Complemented()*/, "PlotSurfaceOffsetOneByMyBspline", 20, 20);

}
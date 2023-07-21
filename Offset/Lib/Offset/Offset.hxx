#pragma once
#include<TopoDS_Shape.hxx>
#include<TopoDS_Edge.hxx>
#include<vector>
#include<gp_Pnt.hxx>
#include<Geom_BSplineCurve.hxx>
#include<Geom_Curve.hxx>
#include<Geom_Surface.hxx>
#include <Geom_Geometry.hxx>
#include<TColgp_Array1OfPnt.hxx>


class Offset
{

private:

	int theNumberOfPoint_;
	std::vector<double> theU_;
	int theNumberOfSection_;
	int theDegree_;


public:


	Offset(int deg);

	std::vector<std::vector<gp_Pnt>> ReadOffsetOne()const;
	std::vector<std::vector<gp_Pnt>> ReadOffsetTwo()const;


	std::vector<Handle(Geom_Curve)> CreatCurveOffsetOneByInterpolate()const;
	std::vector<Handle(Geom_Curve)> CreatCurveOffsetOneBySpline()const;


	std::vector<Handle(Geom_Curve)> CreatCurveOffsetTwoBySpline()const;
	std::vector<Handle(Geom_Curve)> CreatCurveOffsetTwoByInterpolate()const;


	void PlotCurve(std::vector<Handle(Geom_Curve)> curve, std::string filename, int n);


	void PlotSurfaceOffsetOneByGeomFill_BSplineCurves()const;
	void PlotSurfaceOffsetOneByInterpolateCurve()const;
	void PlotSurfaceOffsetOneByOccBsplineCurve()const;
	void PlotSurfaceOffsetOneByMyBspline()const;
	
	
	void PlotSurfaceOffsetOneBySweepApp()const;


	

	void PlotSurfaceOffsetTwoByGeomFill_BSplineCurves()const;
	void PlotSurfaceOffsetTwoByInterpolateCurve()const;




	TopoDS_Shape CreatShapeFromTableOffset1()const;
	TopoDS_Shape CreatShapeFromTableOffset2()const;


	/**************************************************************/
	/****/
	//for separate Bottome and crust->
	std::vector<Handle(Geom_Curve)> CreatCurveTheBottomOftheShip()const;
	std::vector<Handle(Geom_Curve)> CreatCurveTheCrustOftheShip()const;
	


	std::vector<Handle(Geom_Curve)> CreatCurvePathOffsetOneByInterpolateBottomeVec()const;
	



	void PlotCurveTheBottomOftheShip(std::string filename, int n);
	void PlotCurveTheCrustOftheShip(std::string filename, int n);


	void PlotSurfaceTheBottomOftheShipByInterpolateCurve()const;
	void PlotSurfaceTheCrustOftheShipByInterpolateCurve()const;
	void PlotSurfaceOffsetOneTESTBottome()const;


	TopoDS_Shape CreatShapeFromTableOffset1Bottome()const;
	TopoDS_Shape CreatShapeFromTableOffset1Crust()const;

	

	void ShapeOfShip()const;
	
	
	
	
};



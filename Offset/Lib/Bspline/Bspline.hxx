#pragma once
#pragma once
#pragma once
//#include<Shape.hxx>
#include<iostream>
#include<vector>
#include<Geom_Curve.hxx>
#include<TColgp_Array1OfPnt.hxx>
#include<gp_Pnt.hxx>



class Bspline/* : public Shape*/
{

private:

	std::vector<gp_Pnt> thePts_;
	std::vector<double> theU_;
	int theDegree_;

public:

	Bspline() {}

	Bspline(int Deg, std::vector<gp_Pnt> pts);
	void SetDegree(int deg) { theDegree_ = deg; }
	int GetDegree()const { return theDegree_; }
	std::vector<gp_Pnt>GetPts()const { return thePts_; }
	std::vector<double>GetU()const { return theU_; }
	void SetPts(const std::vector<gp_Pnt>& pts) { thePts_ = pts; }
	void SetU(const std::vector<double>& u) { theU_ = u; }

	void CalcTheU();
	void Plot(std::string filename, int n);

	Handle(Geom_Curve) CreateOCCurve()const;

	std::vector<int> CalcMutiplicities()const;
	std::vector<double> CalcKnots() const;


};
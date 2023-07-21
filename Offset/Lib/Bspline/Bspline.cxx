#include"Bspline.hxx"
#include<Geom_BSplineCurve.hxx>
#include<fstream>




Bspline::Bspline(int deg, std::vector<gp_Pnt> pts)
{
	theDegree_ = deg;
	thePts_ = pts;

	CalcTheU();
}


void  Bspline::CalcTheU()
{
	theU_ = CalcKnots();

	return;


	int n = thePts_.size();
	int m = n + theDegree_ + 2;


	for (int j = 0; j < theDegree_ + 1; j++)
	{
		//theU_[j] = 0.0;

		theU_.push_back(0.0);
	}


	for (int j = theDegree_ + 1; j < m - (theDegree_ + 1); j++)
	{
		//theU_[j] = ((1.0) / (m - (2.0 * theDegree_ + 2.0) + 1.0));
		theU_.push_back(theU_[j - 1] + (((1.0) / (m - (2.0 * theDegree_ + 2.0) + 1.0))));
	}

	for (int j = m - (theDegree_ + 1); j < m; j++)
	{
		theU_.push_back(1.0);
	}


}






void Bspline::Plot(std::string filename, int n)
{

	std::fstream My_File(filename, std::ios::out);
	if (!My_File.is_open())
	{
		std::cout << "file not open";
		return;
	}


	My_File << "VARIABLES = X Y Z" << std::endl;
	My_File << "ZONE T = Curve" << std::endl;



	auto curve0 = CreateOCCurve();
	//Curv0.Value()
	double du = 1.0 / n;

	for (int i = 0; i <= n; i++)
	{
		double u = i * du;
		gp_Pnt p = curve0->Value(u);


		My_File << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;

	}

	My_File.close();
}






std::vector<int> Bspline::CalcMutiplicities() const
{
	std::vector<int> mults;

	int sizeOfOther = thePts_.size() - theDegree_ - 1;
	mults.resize(2 + sizeOfOther);

	for (int i = 0; i < mults.size(); i++)
		mults[i] = 1;

	mults[0] = theDegree_ + 1;
	mults[mults.size() - 1] = theDegree_ + 1;

	return mults;
}




std::vector<double> Bspline::CalcKnots() const
{
	std::vector<double> knots;

	knots.resize(2 + thePts_.size() - theDegree_ - 1);
	knots[0] = 0.0;
	knots[knots.size() - 1] = 1.0;

	double dK = 1.0 / (knots.size() - 1.0);

	for (int i = 1; i < knots.size() - 1; i++)
		knots[i] = i * dK;

	return knots;
}



Handle(Geom_Curve) Bspline::CreateOCCurve() const
{
	TColgp_Array1OfPnt Pts;
	Pts.Resize(1, thePts_.size(), false);

	for (int i = 0; i < thePts_.size(); i++)
	{
		double x = thePts_[i].X();
		double y = thePts_[i].Y();
		double z = thePts_[i].Z();

		gp_Pnt p(x, y, z);

		Pts.SetValue(i + 1, p);
	}

	TColStd_Array1OfReal Knot;
	Knot.Resize(1, theU_.size(), false);

	for (int i = 0; i < theU_.size(); i++)
	{
		Knot.SetValue(i + 1, theU_[i]);

	}

	TColStd_Array1OfInteger Mults;

	auto myMults = CalcMutiplicities();
	Mults.Resize(1, myMults.size(), false);

	for (int i = 0; i < myMults.size(); i++)
	{
		Mults.SetValue(i + 1, myMults[i]);
	}

	Handle(Geom_Curve) curve;

	// Get Massage, for ex sqrt(x) 
	try
	{
		curve = new Geom_BSplineCurve(Pts, Knot, Mults, theDegree_, false);
	}
	catch (const Standard_Failure& ex)
	{
		std::cout << ex.GetMessageString() << std::endl;
	}

	return curve;

}


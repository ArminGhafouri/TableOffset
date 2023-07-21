#include"Tools.hxx"
#include<TopExp_Explorer.hxx>
#include<TopoDS.hxx>
#include<TopoDS_Face.hxx>
#include<Geom_BoundedSurface.hxx>
#include<BRep_Tool.hxx>
#include<BRepBuilderAPI_NurbsConvert.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Macro.hxx>
#include <BRep_Builder.hxx>
#include <BRepMesh_PluginEntryType.hxx>
#include <BRepMesh_FactoryError.hxx>
#include <TColStd_MapOfAsciiString.hxx>
#include <BRepMesh_FastDiscret.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TCollection_AsciiString.hxx>
#include <ShapeAnalysis_Edge.hxx>

#include <BRepFilletAPI_MakeFillet.hxx>




void Tools::PlotShape(const TopoDS_Shape& shape, const std::string& filename, int nu, int nv)
{

	std::fstream MyFile(filename + ".plt", std::ios::out | std::ios::app);

	if (!MyFile.is_open())
	{
		std::cout << "file not open" << std::endl;

	}

	//try
	{
		//convert all surfaces to nurbs 
		BRepBuilderAPI_NurbsConvert Shape(shape, false);


		//creat shape



		/// get faces of the shape 
		TopExp_Explorer Ex;

		for (Ex.Init(Shape, TopAbs_FACE); Ex.More(); Ex.Next())
		{
			MyFile << "VARIABLES = X Y Z" << std::endl;
			MyFile << "Zone  I = " << nu + 1 << " ,J = " << nv + 1 << " , Datapacking = Point" << std::endl;

			const TopoDS_Face& Face = TopoDS::Face(Ex.Current());  //-> 
			auto Surface = BRep_Tool::Surface(Face);

			//convert unbunded surface to bounded surface
			//auto BSurface = Handle(Geom_BoundedSurface)::DownCast(Surface);

			double au, bu, av, bv;
			Surface->Bounds(au, bu, av, bv);


			double du = (bu - au) / nu;
			double dv = (bv - av) / nv;


			for (int i = 0; i <= nu; i++)
			{
				double u = i * du + au;


				for (int j = 0; j <= nv; j++)
				{

					double v = j * dv + av;

					gp_Pnt p = Surface->Value(u, v);

					MyFile << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;
				}

			}


		}
	}

	/*catch (const Standard_Failure& ex)
	{
		std::cout << ex.GetMessageString() << std::endl;
	}*/

	MyFile.close();

}



//void Tools::PlotShape2(const TopoDS_Shape& shape, const std::string& filename)
//{
//	std::fstream MyFile(filename + ".plt", std::ios::out | std::ios::app);
//
//	if (!MyFile.is_open())
//	{
//		std::cout << "file not open" << std::endl;
//
//	}
//

//
//
// 
void Tools::PlotShapeTwo(const TopoDS_Shape& shape, const std::string& filename)
{


	std::fstream MyFile(filename + ".plt", std::ios::out /*| std::ios::app*/);

	if (!MyFile.is_open())
	{
		std::cout << "file not open" << std::endl;
	}


	BRepMesh_IncrementalMesh(shape, 0.01);

	//get faces of the shape
	TopExp_Explorer Ex;
	for (Ex.Init(shape, TopAbs_FACE); Ex.More(); Ex.Next())
	{
		const TopoDS_Face& Face = TopoDS::Face(Ex.Current());
		TopLoc_Location loc = Face.Location();
		Handle_Poly_Triangulation Triangul = BRep_Tool::Triangulation(Face, loc);

		if (Triangul.IsNull())
			continue;

		int N = Triangul->NbNodes();
		int E = Triangul->NbTriangles();

		MyFile << "VARIABLES = X Y Z" << std::endl;
		MyFile << "ZONE N =" << N << ", E =" << E << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;



		//GCPnts_TangentialDeflection h;


		//BRepMesh_FastDiscret()
		//FastDiscret(Face);


		for (int i = 1; i <= Triangul->NbNodes(); i++)
		{

			double x = Triangul->Node(i).X();
			double y = Triangul->Node(i).Y();
			double z = Triangul->Node(i).Z();

			MyFile << x << "  " << y << "  " << z << std::endl;

		}

		for (int i = 1; i <= E; i++)
		{
			auto t = Triangul->Triangle(i);
			int N1, N2, N3;
			t.Get(N1, N2, N3);
			MyFile << N1 << "  " << N2 << "  " << N3 << std::endl;

		}

	}

	MyFile.close();  //////*****>>>>>>>>>>>

}




TopoDS_Shape Tools::JoinOCShape(std::vector<TopoDS_Face> face)
{

	TopoDS_Compound FinallShape;
	BRep_Builder aBuilder;

	aBuilder.MakeCompound(FinallShape);

	for (int i = 0; i < face.size(); i++)
	{
		aBuilder.Add(FinallShape, face[i]);
	}

	return FinallShape;

}



//TopoDS_Shape Tools::MakeFillet(TopoDS_Shape shape, double FilletSize)
//{
//
//	BRepFilletAPI_MakeFillet mkFillet(shape);
//	TopExp_Explorer anEdgeExplorer(shape, TopAbs_EDGE);
//
//	{
//		TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
//
//
//		mkFillet.Add(FilletSize, anEdge);
//		anEdgeExplorer.Next();
//
//
//	}
//
//	TopoDS_Shape shape0 = mkFillet.Shape();
//
//	return shape0;
//}





//TopoDS_Face Tools::JoinOCEdge(std::vector<TopoDS_Vertex> vertex)
//{
//
//	TopoDS_Face FinallFace;
//	BRep_Builder aBuilder;
//
//	aBuilder.MakeFace(FinallFace);
//
//	for (int i = 0; i < vertex.size(); i++)
//	{
//		aBuilder.Add(FinallFace, vertex[i]);
//
//	}
//
//
//	return FinallFace;
//
//}




std::vector <int> Tools::FindeIntegerOfEdges(TopoDS_Shape shape, std::vector<TopoDS_Edge>edgs)
{
	std::vector<int> Numbers;

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(shape, TopAbs_EDGE, mymap);

	for (int i = 0; i < edgs.size(); i++)
	{



		int number = mymap.FindIndex(edgs[i]);
		Numbers.push_back(number);

		std::cout << "edgs" << i << "{" << number << "}" << std::endl;

	}

	return Numbers;

}




std::vector<TopoDS_Edge> Tools::FindEdges(TopoDS_Shape shape, std::vector<int> index)
{

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(shape, TopAbs_EDGE, mymap);
	std::vector<TopoDS_Edge> edgs;


	for (int i = 0; i < index.size(); i++)
	{
		TopoDS_Edge edge = TopoDS::Edge(mymap(index[i]));
		edgs.push_back(edge);
	}

	return edgs;
}



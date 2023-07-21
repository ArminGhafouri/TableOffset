#pragma once
#pragma once
#include<TopoDS_Shape.hxx>
#include<TopoDS_Edge.hxx>
#include<iostream>
#include<vector>
#ifndef Tools_HEADER
#define Tools_HEADER
#include <Poly_Triangulation.hxx>


class TopoDS_Face;

class Tools
{

private:


public:

	static void PlotShape(const TopoDS_Shape& shape, const std::string& filename, int nu, int nv);


	static void PlotShapeTwo(const TopoDS_Shape& shape, const std::string& filename);
	static TopoDS_Shape JoinOCShape(std::vector<TopoDS_Face> face);

	/*static TopoDS_Face JoinOCEdge(std::vector<TopoDS_Vertex> vertex);*/
	//static TopoDS_Shape MakeFillet(TopoDS_Shape shape, double FilletSize);


	static std::vector<int> FindeIntegerOfEdges(TopoDS_Shape shape, std::vector<TopoDS_Edge> edgs);
	static std::vector<TopoDS_Edge> FindEdges(TopoDS_Shape shape, std::vector<int> index);


};

#endif
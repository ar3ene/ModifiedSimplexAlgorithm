#ifndef _NELDER_MEAD_ALGORITHM_H_
#define _NELDER_MEAD_ALGORITHM_H_

#include <iostram>
#include <vector>
#include "sys_types.h"

#pragma once

class CNelderMeadAlgorithm
{
public:
	CNelderMeadAlgorithm(int num_params, int num_vertices, int trace_path);
	~CNelderMeadAlgorithm();

	void Init(const std::vector<std::vector<double>>& initial_simplex);

	void Curvefitting0(FUNCTION function, std::vector<double>& optimal,
		int num_iterations, double tolerance, int num_params,
		double reflection_coeff, double expansion_coeff, double contraction_coeff);
	
	void Curvefitting1(FUNCTION function, std::vector<double>& optimal,
		int num_iterations, double tolerance, int num_params,
		double reflection_coeff, double expansion_coeff, double contraction_coeff);
	
	void Curvefitting2(FUNCTION function, std::vector<double>& optimal,
		int num_iterations, double tolerance, int num_params,
		double reflection_coeff, double expansion_coeff, double contraction_coeff);
	
	void Curvefitting3(FUNCTION function, std::vector<double>& optimal,
		int num_iterations, double tolerance, int num_params,
		double reflection_coeff, double expansion_coeff, double contraction_coeff);

private:
	void Reflect(Point& centroid, Point& p, Point& p_reflected, 
		double reflection_coeff);

	void Expand(std::vector<double>& centroid, std::vector<double>& p_expanded, 
		std::vector<double>& p_reflected, double expansion_coeff);

	void ContractOneDim(std::vector<std::vector<double>>& vertexes_ensemble, 
		std::vector<double>& centroid, std::vector<double>& p_contracted, 
		int max_idx, double contraction_coeff);

	void ContractMultiDim(std::vector<std::vector<double>>& vertexes_ensemble, 
		int min_idx);

	void ExpandMultiDim(int min_idx, double expansion_coeff);

	void ComputeVertexValues(std::vector<std::vector<double>>& vertexes_ensemble, 
		std::vector<double>& f_values, FUNCTION Function);

	void ComputeCentroid(std::vector<std::vector<double>>& vertexes_ensemble,
		std::vector<double>& centroid, int num_params, int max_idx);

	void MiddlePoint(const Point& a, const Point& b, Point& mid);

	void RelativeAdd(const Point& o, const Point& a, const Point& b, Point& s);

	void RelativeSub(const Point& o, const Point& a, const Point& b, Point& s);

	int FindMinIndex(std::vector<double>& values);

	int FindMaxIndex(std::vector<double>& values);

	bool IsConverged(double tolerance);

	void TraceNewMin(Point& p, int num_params, double fmin, int idx, const char* info);

	std::vector<std::vector<double>> m_simplex;

	int m_num_params;
	int m_num_vertexes;
	int m_trace_path;
};

#endif /// ! _NELDER_MEAD_ALGORITHM_H_ 
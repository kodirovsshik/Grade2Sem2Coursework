
#include <vector>
#include <iostream>
#include <unordered_set>

#include <cassert>


using weights_t = std::vector<double>;

using upper_limit_t = double;
struct constraint_t
{
	weights_t weights;
	upper_limit_t limit;
};

using values_t = std::vector<double>;


double evaluate(const weights_t& coeffs, const values_t& values)
{
	assert(coeffs.size() == values.size());

	double sum = 0;
	for (size_t i = 0; i < coeffs.size(); ++i)
		sum += coeffs[i] * values[i];
	return sum;
}

values_t _maximizing_simplex_method_phase2(const std::vector<constraint_t>& constraints, const weights_t& objective)
{
	//vars:
	//[0; N-1] = in/out vars
	//[N; N+M-1] = dummy vars
	//[N+M] = objective

	//tableau:
	//         (IO vars) (dummy vars) (objective) | RHS
	// eq1
	// ...
	// eqM

	static constexpr double error_tolerance = 1e-8;

	using row_t = std::vector<double>;
	using matrix_t = std::vector<row_t>;

	const size_t N = objective.size(); //IO vars
	const size_t M = constraints.size(); //dummy vars

	for (const auto& constraint : constraints)
		assert(constraint.weights.size() == N);

	const size_t width = N + M + 2;
	const size_t vars_count = N + M;
	const size_t height = M + 1;



	std::vector<size_t> row_associated_vars;
	row_associated_vars.resize(M);
	for (size_t i = 0; i < M; ++i)
		row_associated_vars[i] = i + N;



	matrix_t tableau;
	tableau.resize(height);
	
	for (size_t i = 0; i < M; ++i)
	{
		tableau[i].resize(width);
		for (size_t j = 0; j < N; ++j)
			tableau[i][j] = constraints[i].weights[j];
		for (size_t j = N; j < N + M; ++j)
			tableau[i][j] = double(i == (j - N));
		tableau[i][N + M] = 0;
		tableau[i][N + M + 1] = constraints[i].limit;
	}

	tableau[M].resize(width);
	for (size_t j = 0; j < N; ++j)
		tableau[M][j] = -objective[j];
	for (size_t j = N; j < N + M; ++j)
		tableau[M][j] = 0;
	tableau[M][N + M] = 1;
	tableau[M][N + M + 1] = 0;
	


	auto get_simplex_ratio = [&]
	(size_t i, size_t j) -> double
	{
		if (tableau[i][j] <= 0)
			return INFINITY;
		return  tableau[i][N + M + 1] / tableau[i][j];
	};


	auto try_rowwise_transformation = [&]
	(size_t min_column) -> bool
	{
		size_t min_row = 0;
		double min_ratio = get_simplex_ratio(0, min_column);
		for (size_t i = 1; i < M; ++i)
		{
			double current_ratio = get_simplex_ratio(i, min_column);
			if (current_ratio < min_ratio)
			{
				min_ratio = current_ratio;
				min_row = i;
			}
		}

		if (min_ratio == INFINITY)
			return false;

		const size_t old_basis = row_associated_vars[min_row];
		if (old_basis == min_column)
			return false;

		const double intersection_elem = tableau[min_row][min_column];
		for (size_t j = 0; j < width; ++j)
			tableau[min_row][j] /= intersection_elem;

		for (size_t i = 0; i <= M; ++i)
		{
			if (i == min_row)
				continue;

			const double row_factor = tableau[i][min_column];
			for (size_t j = 0; j < width; ++j)
				tableau[i][j] -= tableau[min_row][j] * row_factor;
		}

		row_associated_vars[min_row] = min_column;
		return true;
	};

	size_t iters = vars_count * vars_count + 1;
	while (iters --> 0)
	{
		auto& objective_row = tableau[M];
		size_t min_column = 0;
		for (size_t i = 1; i < vars_count; ++i)
		{
			if (objective_row[i] < objective_row[min_column])
				min_column = i;
		}

		if (objective_row[min_column] > error_tolerance)
			break;

		if (objective_row[min_column] >= -error_tolerance)
		{
			break;
		}
			//got zero as best choice, try finding at least something
		//	static std::vector<size_t> idx;
		//	idx.clear();
		//	idx.resize(vars_count);
		//	for (size_t i = 0; i < vars_count; ++i)
		//		idx[i] = i;

		//	bool ok = false;
		//	while (!idx.empty())
		//	{
		//		const size_t n = rand() % idx.size();
		//		const size_t j = idx[n];
		//		idx.erase(idx.begin() + n);

		//		if (objective_row[j] > error_tolerance)
		//			continue;
		//		if (true == (ok = try_rowwise_transformation(j)))
		//			break;
		//	}
		//	if (!ok)
		//		break;
		//}
		//else
		//{
		//	if (!try_rowwise_transformation(min_column))
		//		break;
		//}

	}



	values_t result(N, 0);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < M; ++j)
		{
			if (row_associated_vars[j] == i)
			{
				result[i] = tableau[j][N + M + 1];
				break;
			}
		}
	}

	for (const auto& constraint : constraints)
	{
		const double a = evaluate(constraint.weights, result);
		const double b = constraint.limit;
		const double max_err = fabs(a + b) / 2 * error_tolerance;
		const double diff = a - b;
		if (diff > max_err)
			return {};
	}

	return result;
}


values_t maximizing_simplex_method(std::vector<constraint_t> constraints, const weights_t& objective)
{

}

values_t minimizing_simplex_method(std::vector<constraint_t> constraints, weights_t objective)
{
	//for (auto& x : objective)
	//	x = -x;

	//const size_t n = objective.size();
	//weights_t tmp_obj(n * 2);
	//for (size_t i = 0; i < n; ++i)
	//	tmp_obj[i + n] = 1;

	//return _maximizing_simplex_method_phase2(constraints, objective);
}
values_t minimizing_integer_simplex_method(const std::vector<constraint_t>& constraints, const weights_t& objective)
{
	return {};
}




int main()
{
	const std::vector<constraint_t> constraints = {
		{{1, 0}, 3},
		{{-1, 1}, 1},
		{{-1, 0}, -1},
	};

	const weights_t objective{-2, 1};

	auto solution = maximizing_simplex_method(constraints, objective);
	

	return 0;
}

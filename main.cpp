
#include <vector>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <ranges>
#include <source_location>

#include <clocale>

#include <conio.h>



using weights_t = std::vector<double>;

using upper_limit_t = double;
struct constraint_t
{
	weights_t weights;
	upper_limit_t limit = NAN;
	int direction = 1;
};
using constraints_t = std::vector<constraint_t>;

using values_t = std::vector<double>;



#define xassert(cond) _xassert((cond), #cond)
void _xassert(bool ok, const char* cond, std::source_location pos = std::source_location::current())
{
	if (!ok)
	{
		std::cerr << "\aFrom " << pos.file_name() << ':' << pos.line() << ":\n";
		std::cerr << " In function " << pos.function_name() << ": ASSERTION FAILED:\n";
		std::cerr << "  " << cond << std::endl;
		system("pause");
	}
}



void free(values_t& x)
{
	values_t empty;
	x.swap(empty);
}


double evaluate(const weights_t& coeffs, const values_t& values)
{
	xassert(coeffs.size() == values.size());

	double sum = 0;
	for (size_t i = 0; i < coeffs.size(); ++i)
		sum += coeffs[i] * values[i];
	return sum;
}


static constexpr double error_tolerance = 1e-4;
double get_error_tolerance(double around)
{
	around = fabs(around);
	if (around < error_tolerance)
		return error_tolerance;
	return error_tolerance * around;
}

values_t minimizing_simplex_method(const constraints_t& constraints, const weights_t& objective)
{
	//vars:
	//[0; N-1] = vars
	//[N] = objective
	//[N+1; N+M] = dummy vars

	//tableau:
	//            (vars) (obj) (dummy vars) [obj2] (artificial vars) | RHS
	//(ineqs)
	//(obj)
	//[obj2]

	//where obj2 is an artificial objective for the purpose of canonnicalizing


	using row_t = std::vector<double>;
	using matrix_t = std::vector<row_t>;

	size_t n_vars = objective.size();
	size_t n_constraints = constraints.size();
	size_t n_equational_vars = n_constraints;
	size_t n_artificial_vars = n_constraints;

	for (const auto& constraint : constraints)
		xassert(constraint.weights.size() == n_vars);

	size_t total_vars_count = n_vars + n_equational_vars + 1 + n_artificial_vars;
	size_t width = total_vars_count + 2;
	size_t height = n_constraints + 2;
	size_t objective_column = width - 1;


	std::vector<size_t> row_associated_vars;
	row_associated_vars.resize(height, -1);



	matrix_t tableau;
	tableau.resize(height);



	auto get_simplex_ratio = [&]
	(size_t i, size_t j) -> double
	{
		if (tableau[i][j] <= 0)
			return INFINITY;
		return  tableau[i][width - 1] / tableau[i][j];
	};

	auto add_row = [&]
	(size_t to, size_t which, double factor)
	{
		for (size_t j = 0; j < width; ++j)
			tableau[to][j] += tableau[which][j] * factor;
	};

	auto try_rowwise_transformation = [&]
	(size_t pivot_column) -> bool
	{
		size_t min_row = 0;
		double min_ratio = get_simplex_ratio(min_row, pivot_column);
		for (size_t i = 1; i < n_constraints; ++i)
		{
			double current_ratio = get_simplex_ratio(i, pivot_column);
			if (current_ratio < min_ratio)
			{
				min_ratio = current_ratio;
				min_row = i;
			}
		}

		if (min_ratio == INFINITY)
			return false;

		const size_t old_basis = row_associated_vars[min_row];
		if (old_basis == pivot_column)
			return false;

		const double intersection_elem = tableau[min_row][pivot_column];
		for (size_t j = 0; j < width; ++j)
			tableau[min_row][j] /= intersection_elem;

		for (size_t i = 0; i < height; ++i)
		{
			if (i == min_row)
				continue;

			add_row(i, min_row, -tableau[i][pivot_column]);
		}

		row_associated_vars[min_row] = pivot_column;
		return true;
	};

	auto pick_pivoting_column = [&](size_t start) -> size_t
	{
		if (start >= total_vars_count - 1)
			return -1;

		size_t pivot_column = start;
		auto& objective_row = tableau[height - 1];
		for (size_t i = start + 1; i < total_vars_count; ++i)
		{
			if (i == objective_column)
				continue;
			if (objective_row[i] > objective_row[pivot_column])
				pivot_column = i;
		}
		return objective_row[pivot_column] >= 0 ? pivot_column : -1;
	};

	auto solve_current_objective = [&]
	{
		auto& vars = row_associated_vars;
		size_t iters = total_vars_count * total_vars_count;
		size_t starting_column = 0;
		while (iters-- > 0)
		{
			auto& objective_row = tableau[height - 1];
			size_t pivot_column = pick_pivoting_column(starting_column);

			if (pivot_column == -1)
				break;

			if (objective_row[pivot_column] <= error_tolerance)
				break;

			if (!try_rowwise_transformation(pivot_column))
			{
				++iters;
				starting_column = pivot_column + 1;
				continue;
			}

			starting_column = 0;
		}
	};



	//Constraints data
	for (size_t i = 0; i < n_constraints; ++i)
	{
		tableau[i].resize(width);

		for (size_t j = 0; j < n_vars; ++j)
			tableau[i][j] = constraints[i].weights[j];
		tableau[i][width - 1] = constraints[i].limit;

		tableau[i][n_vars + 1 + i] = constraints[i].direction;

		if (constraints[i].limit < 0)
			add_row(i, i, -2); //negate the row

		tableau[i][n_vars + 1 + n_equational_vars + 1 + i] = 1;
	}


	//Original objective row
	tableau[height - 2].resize(width);
	for (size_t j = 0; j < n_vars; ++j)
		tableau[height - 2][j] = -objective[j];
	tableau[height - 2][n_vars] = 1;


	//Artificial objective row
	tableau[height - 1].resize(width);
	for (size_t j = 0; j < n_equational_vars; ++j)
		tableau[height - 1][n_vars + 1 + n_equational_vars + 1 + j] = -1;
	tableau[height - 1][n_vars + 1 + n_equational_vars] = 1;

	auto artificial_objective = [&]
	{
		return fabs(tableau[height - 1][width - 1]);
	};


	for (size_t i = 0; i < n_constraints; ++i)
		add_row(height - 1, i, 1);


	{ int _ = 0; }

	//Phase 1: cannonicalizing
	for (size_t i = 0; i < n_vars + 1 + n_equational_vars; ++i)
		try_rowwise_transformation(i);

	if (artificial_objective() > error_tolerance)
	{
		solve_current_objective();
		if (artificial_objective() > error_tolerance)
			return {};
	}


	--height;
	width -= n_artificial_vars + 1;
	total_vars_count = n_vars + 1 + n_equational_vars;
	objective_column = n_vars;

	tableau.pop_back();
	for (auto& row : tableau)
	{
		row[n_vars + 1 + n_equational_vars] = row.back();
		row.resize(width);
	}



	//Phase 2: solving objective
	solve_current_objective();



	values_t result(n_vars, 0);
	for (size_t i = 0; i < n_vars; ++i)
	{
		for (size_t j = 0; j < n_constraints; ++j)
		{
			if (row_associated_vars[j] == i)
			{
				result[i] = tableau[j][width - 1];
				break;
			}
		}
	}

	for (const auto& constraint : constraints)
	{
		if (constraint.direction == 0)
		{
			const double a = evaluate(constraint.weights, result);
			const double b = constraint.limit;
			const double max_err = fabs(a + b) / 2 * error_tolerance;
			if (fabs(a - b) > max_err)
				return {};
		}
		else
		{
			const double a = evaluate(constraint.weights, result);
			const double b = constraint.limit;
			const double max_err = fabs(a + b) / 2 * error_tolerance;
			const double diff = a - b;
			const double k = constraint.direction;
			if (diff * k > max_err)
				return {};
		}
	}

	result.push_back(tableau[height - 1][width - 1]);
	return result;
}

values_t maximizing_simplex_method(const constraints_t& constraints, weights_t objective)
{
	for (auto& x : objective)
		x = -x;

	auto result = minimizing_simplex_method(constraints, objective);
	if (result.size())
		result.back() *= -1;

	return result;
}

values_t minimizing_integer_simplex_method(constraints_t& constraints, const weights_t& objective)
{
	auto solution = minimizing_simplex_method(constraints, objective);
	if (solution.empty())
		return {};

	const double optimal_value = solution.back();
	solution.pop_back();

	auto integer_divergence = []
	(double x) -> double
	{
		return fabs(x - round(x));
	};

	const size_t most_divergent_index = std::ranges::max_element(solution, {}, integer_divergence) - solution.begin();
	const double most_divergent = solution[most_divergent_index];

	if (integer_divergence(most_divergent) < get_error_tolerance(most_divergent))
	{
		solution.push_back(optimal_value);
		return solution;
	}

	free(solution);

	constraint_t new_constraint;
	new_constraint.weights.resize(objective.size());

	//x[i] <= floor(X[i]);
	new_constraint.weights[most_divergent_index] = 1;
	new_constraint.direction = 1;
	new_constraint.limit = floor(most_divergent);

	constraints.push_back(std::move(new_constraint));
	auto solution_left = minimizing_integer_simplex_method(constraints, objective);
	

	new_constraint = std::move(constraints.back());
	new_constraint.direction = -1;
	new_constraint.limit += 1;
	constraints.back() = std::move(new_constraint);

	auto solution_right = minimizing_integer_simplex_method(constraints, objective);

	constraints.pop_back();

	if (solution_left.empty())
		return solution_right;
	if (solution_right.empty())
		return solution_left;

	if (solution_right.back() < solution_left.back())
		return solution_right;
	else
		return solution_left;
}
values_t minimizing_integer_simplex_method(const constraints_t& constraints, const weights_t& objective)
{
	constraints_t copy = constraints;
	return minimizing_integer_simplex_method(copy, objective);
}


struct task_t
{
	double price_new = 400;
	double price_old = 800;

	double hours_new = 50;
	double hours_old = 100;

	double gamma = 0.9;

	double initial_workers = 60;

	double min_hours[6] = { 6500, 6600, 6600, 6500, 6400, 6400 };



	double gamma_sum(size_t n)
	{
		if (gamma == 1)
			return (double)n;
		else
			return (1 - pow(gamma, n + 1)) / (1 - gamma);
	};

	void solve()
	{
		std::vector<constraint_t> constraints;

		for (size_t i = 1; i <= 6; ++i)
		{
			constraint_t c;
			c.direction = -1;
			c.limit = min_hours[i - 1] - initial_workers * pow(gamma, i - 1) * hours_old;

			c.weights.resize(6, 0);
			for (size_t j = 1; j <= i - 1; ++j)
				c.weights[j - 1] = pow(gamma, i - 1 - j) * hours_old;
			c.weights[i - 1] = hours_new;

			constraints.push_back(std::move(c));
		}



		weights_t objective(6, 0);

		for (size_t i = 1; i <= 6; ++i)
		{
			objective[i - 1] = price_new + price_old * gamma_sum(5 - i);
		}

		auto solution = minimizing_integer_simplex_method(constraints, objective);
		if (solution.size())
			solution.back() += initial_workers * price_old * gamma_sum(5);

		return solution;
	}
};

task_t task{};

void editor()
{
	while (true)
	{
		std::wcout << L"\n Текущие данные:\n";

		std::wcout << L" [1] З/П стажёра: " << task.price_new << '\n';
		std::wcout << L" [2] З/П работника: " << task.price_old << '\n';
		std::wcout << L" [3] Кол-во часов стажёра в месяц: " << task.hours_new << '\n';
		std::wcout << L" [4] Кол-во часов работника в месяц: " << task.hours_old << '\n';
		std::wcout << L" [5] Максимальный процент увольняемости: " << 100 * (1 - task.gamma) << '\n';
		std::wcout << L" [6] Изначальное количество работников: " << task.min_hours << '\n';
		std::wcout << L" [7] Необходимое число рабочих часов в месяц:\n  ";

		for (auto h : task.min_hours)
			std::wcout << h << ' ';

		std::wcout << "\n\n";
		std::wcout << L" Нажмите на номер переменной или массива для изменения значения\n";
		std::wcout << L" Нажмите на 0 чтобы вернуться в меню\n";
		std::wcout << L" Нажмите на Delete чтобы восстановить исходные значения\n";
		std::wcout << L"\n >>> ";

		char c = _getch();
		if (c == -1)
			task = task_t();
		else if (c == '0')
			return;
		else
		{
			if (c < '1' || c > '7')
				continue;

			std::wcout << c << "\n\n";
			c -= '0';
		
			double* p_var = nullptr;
			double dummy = 0;
			if (c == 7)
			{
				int x = 0;
				while (x < 1 || x > 7)
				{
					std::wcout << " Введите номер месяца (1 - Январь, 6 - Июнь)\n";
					std::wcout << L" >>> ";
					std::cin >> x;
				}
				p_var = &task.min_hours[x - 1];
			}
			else if (c == 5)
			{

			}
		}
	}
}

void about()
{
	std::wcout << L"\n\
 Данная программа была разработана студентом группы А01ИСТ2 МГЭИ А. Д. Сахарова \
 БГУ Егоровым Станиславом в рамках курсовой работы по дисциплине \"Исследование \n\
 операций\"\n\n\
 Вариант: 18\n\
 Тема: \"Задача линейного программирования\"\n\n\
 Нажмите любую клавишу для продолжения";
	(void)_getch();
}

void menu()
{
	static const wchar_t menu_data[] = 
LR"(
 Меню:
 1 - Решить задачу
 2 - Напечатать/редактировать входные данные
 3 - Об авторе
 4 - Выход
 >>> )";
	while (true)
	{
		system("cls");
		std::wcout << menu_data;

		int choice;
		std::cin >> choice;
		
		system("cls");
		switch (choice)
		{
		case 1:
			task.solve();
			break;

		case 2:
			//editor();
			break;

		case 3:
			about();
			break;

		case 4:
			return;
		}
 	}
}

int main()
{
	xassert(setlocale(LC_ALL, ".65001"));

	menu();

	return 0;
}
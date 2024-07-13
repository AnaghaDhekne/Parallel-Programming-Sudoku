#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <mutex>

using namespace std;

using Position = pair<int, int>;

// Abstract base class for Sudoku solvers
class SudokuSolver {
protected:
	vector<vector<int>> board;  // Sudoku board
	int boardSize;              // Size of the board (e.g., 9x9)
	int boxSize;                // Size of each box (e.g., 3x3)
	bool solved;                // Flag to indicate if the puzzle is solved
	int emptyCells;             // Number of empty cells in the puzzle
	int emptyCellsValue;        // Value used to represent empty cells
	int minValue;               // Minimum value allowed in the puzzle
	int maxValue;               // Maximum value allowed in the puzzle

public:
	// Constructor to initialize the Sudoku solver with input from a file
	SudokuSolver(const string& filename) {
		// Open input file
		ifstream inputFile(filename);
		if (!inputFile) {
			cerr << "Error opening file " << filename << endl;
			exit(1);
		}

		// Read board size and initialize board
		inputFile >> boardSize;
		boxSize = sqrt(boardSize);
		board.resize(boardSize, vector<int>(boardSize));

		// Read board values and count empty cells
		int numEmptyCells = 0;
		for (int row = 0; row < boardSize; ++row) {
			for (int col = 0; col < boardSize; ++col) {
				int value;
				inputFile >> value;
				board[row][col] = value;
				numEmptyCells += (value == 0);
			}
		}
		emptyCells = numEmptyCells;
		emptyCellsValue = 0;
		minValue = 1;
		maxValue = boardSize;
		solved = false;

		inputFile.close();
	}

	// Pure virtual function to solve the Sudoku puzzle (implemented by derived classes)
	virtual bool solve() = 0;

	void printSolution() const {
		// Print the solved Sudoku board
		for (int i = 0; i < boardSize; ++i) {
			if (i % boxSize == 0 && i != 0) {
				cout << string(boardSize * 3 + boxSize - 1, '-') << endl;
			}

			for (int j = 0; j < boardSize; ++j) {
				if (j % boxSize == 0 && j != 0) {
					cout << "| ";
				}

				cout << setw(2) << board[i][j] << " ";
			}
			cout << endl;
		}
	}

	void writeOutput(const string& filename, long long duration) const {
		// Write the solved Sudoku board and time taken to an output file
		ofstream outputFile(filename);

		int digit = log10(boardSize) + 1;

		for (int r = 0; r < boardSize; ++r) {
			for (int c = 0; c < boardSize; ++c) {
				outputFile << setw(digit) << board[r][c];

				if (c != boardSize - 1) {
					outputFile << " ";
				}

				if (c % boxSize == (boxSize - 1) && c != boardSize - 1) {
					outputFile << "  ";
				}
			}

			if (r != boardSize - 1) {
				outputFile << "\n";
				if (r % boxSize == (boxSize - 1)) {
					outputFile << "\n";
				}
			}
		}

		outputFile << "\n\n";
		outputFile << "Time taken: " << duration << " microseconds";

		outputFile.close();
	}
};

// Derived class for solving Sudoku using brute force approach
class SudokuSolverBruteForce : public SudokuSolver {
private:
	bool isValidRow(int num, Position pos) const
	{
		// Check if the number is valid in the given row
		for (int i = 0; i < boardSize; ++i)
		{
			if ((i != pos.second) && (board[pos.first][i] == num)) { return false; }
		}

		return true;
	}

	bool isValidColumn(int num, Position pos) const
	{
		// Check if the number is valid in the given column
		for (int i = 0; i < boardSize; ++i)
		{
			if ((i != pos.first) && (board[i][pos.second] == num)) { return false; }
		}

		return true;
	}

	bool isValidBox(int num, Position pos) const
	{
		// Check if the number is valid in the given box
		int box_x = pos.first - pos.first % boxSize;
		int box_y = pos.second - pos.second % boxSize;

		for (int i = box_x; i < box_x + boxSize; ++i)
		{
			for (int j = box_y; j < box_y + boxSize; ++j)
			{
				if ((i != pos.first || j != pos.second) && (board[i][j] == num)) { return false; }
			}
		}

		return true;
	}

	bool isValid(int num, Position pos) const
	{
		// Check if the number is valid in the given position
		return isValidRow(num, pos) && isValidColumn(num, pos) && isValidBox(num, pos);
	}

	// Set a number at a given position on the board
	void set_board_data(int row, int col, int num) { board[row][col] = num; }

public:
	// Constructor to initialize the Sudoku solver with input from a file
	SudokuSolverBruteForce(const string& filename) : SudokuSolver(filename) {}

	// Solve the Sudoku puzzle using brute force approach
	bool solve() override {
		for (int row = 0; row < boardSize; ++row) {
			for (int col = 0; col < boardSize; ++col) {
				if (board[row][col] == emptyCellsValue) {
					for (int num = minValue; num <= maxValue; ++num) {
						Position pos = make_pair(row, col);

						if (isValid(num, pos)) {
							set_board_data(row, col, num);

							if (solve()) {
								solved = true;
								return true;
							}
							else {
								set_board_data(row, col, emptyCellsValue);
							}
						}
					}
					return false;
				}
			}
		}

		solved = true;
		return true;
	}
};

// Derived class for solving Sudoku using backtracking approach
class SudokuSolverBackTracking : public SudokuSolver {
private:
	bool isValidRow(int num, Position pos) const
	{
		// Check if the number is valid in the given row
		for (int i = 0; i < boardSize; ++i)
		{
			if ((i != pos.second) && (board[pos.first][i] == num)) { return false; }
		}

		return true;
	}

	bool isValidColumn(int num, Position pos) const
	{
		// Check if the number is valid in the given column
		for (int i = 0; i < boardSize; ++i)
		{
			if ((i != pos.first) && (board[i][pos.second] == num)) { return false; }
		}

		return true;
	}

	bool isValidBox(int num, Position pos) const
	{
		// Check if the number is valid in the given box
		int box_x = pos.first - pos.first % boxSize;
		int box_y = pos.second - pos.second % boxSize;

		for (int i = box_x; i < box_x + boxSize; ++i)
		{
			for (int j = box_y; j < box_y + boxSize; ++j)
			{
				if ((i != pos.first || j != pos.second) && (board[i][j] == num)) { return false; }
			}
		}

		return true;
	}

	bool isValid(int num, Position pos) const
	{
		// Check if the number is valid in the given position
		return isValidRow(num, pos) && isValidColumn(num, pos) && isValidBox(num, pos);
	}

	// Set a number at a given position on the board
	void set_board_data(int row, int col, int num) { board[row][col] = num; }

	bool checkIfAllFilled() const
	{
		// Check if all cells in the board are filled
		for (int i = 0; i < boardSize; ++i)
		{
			for (int j = 0; j < boardSize; ++j)
			{
				if (board[i][j] == emptyCellsValue)
				{
					return false;
				}
			}
		}
		return true;
	}

	Position find_empty () const
	{
		// Find the next empty cell in the board
		Position empty_cell;
		bool stop = false;

		for (int i = 0; i < boardSize; ++i)
		{
			for (int j = 0; j < boardSize; ++j)
			{
				if (board[i][j] == emptyCellsValue)
				{
					empty_cell = make_pair(i, j);
					stop = true;
					break;
				}
			}
			if (stop) { break; }
		}

		return empty_cell;  // (row, col)
	}

	vector<int> getNumbersInRow(int indexOfRows) const
	{
		// Get all the numbers present in the given row
		vector<int> numbersInRow;

		for (int col = 0; col < boardSize; ++col)
		{
			int num = board[indexOfRows][col];
			if (num == emptyCellsValue) continue;
			numbersInRow.push_back(num);
		}

		return numbersInRow;
	}

	vector<int> getNumbersInCol(int indexOfColumns) const
	{
		// Get all the numbers present in the given column
		vector<int> numbersInCol;

		for (int row = 0; row < boardSize; ++row)
		{
			int num = board[row][indexOfColumns];
			if (num == emptyCellsValue) continue;
			numbersInCol.push_back(num);
		}

		return numbersInCol;
	}

	bool isUnique(int num, Position pos) const
	{
		// Check if the number is unique in its row, column, and box
		int local_row = pos.first % boxSize;
		int local_col = pos.second % boxSize;

		int box_x = floor(pos.first / boxSize);
		int box_y = floor(pos.second / boxSize);

		for (int i = ((local_row == 0) ? 1 : 0); i < boxSize; ++i)
		{
			if (i == local_row) { continue; }
			vector<int> numbersInRow = getNumbersInRow(box_x * boxSize + i);
			if (find(numbersInRow.begin(), numbersInRow.end(), num) == numbersInRow.end())
			{
				return false;
			}
		}

		for (int j = ((local_col == 0) ? 1 : 0); j < boxSize; ++j)
		{
			if (j == local_col) { continue; }
			vector<int> numbersInCol = getNumbersInCol(box_y * boxSize + j);
			if (find(numbersInCol.begin(), numbersInCol.end(), num) == numbersInCol.end())
			{
				return false;
			}
		}

		return true;
	}

	public:
		// Constructor to initialize the Sudoku solver with input from a file
		SudokuSolverBackTracking(const string& filename) : SudokuSolver(filename) {}

		bool solve() override {
			// Solve the Sudoku using backtracking approach
			if (solved) {
				return true;
			}

			if (checkIfAllFilled()) {
				solved = true;
				return true;
			}

			Position emptyPos = find_empty();
			int row = emptyPos.first;
			int col = emptyPos.second;

			for (int num = minValue; num <= maxValue; ++num) {
				if (isValid(num, emptyPos)) {
					set_board_data(row, col, num);

					if (isUnique(num, emptyPos)) {
						num = maxValue + 1;
					}

					if (solve()) {
						solved = true;
						return true;
					}
					else {
						set_board_data(row, col, emptyCellsValue);
					}
				}
			}

			solved = false;
			return false;
		}
};

//class SudokuSolverParallel : public SudokuSolver {
//private:
//	mutex mtx;
//
//	bool is_valid(int row, int col, int num) {
//		int box_row = (row / boxSize) * boxSize;
//		int box_col = (col / boxSize) * boxSize;
//
//		for (int i = 0; i < boardSize; i++) {
//			unique_lock<mutex> lock(mtx);
//			if (board[row][i] == num || board[i][col] == num ||
//				board[box_row + i / boxSize][box_col + i % boxSize] == num)
//				return false;
//		}
//		return true;
//	}
//
//	pair<int, int> find_empty() {
//		for (int i = 0; i < boardSize; i++) {
//			for (int j = 0; j < boardSize; j++) {
//				unique_lock<mutex> lock(mtx);
//				if (board[i][j] == 0)
//					return { i, j };
//			}
//		}
//		return { -1, -1 };
//	}
//
//	bool solve_sudoku(int& empty_count) {
//		auto find = find_empty();
//		if (find.first == -1)
//			return true;
//
//		int row = find.first;
//		int col = find.second;
//
//		for (int i = 1; i <= boardSize; i++) {
//			if (is_valid(row, col, i)) {
//				unique_lock<mutex> lock(mtx);
//				board[row][col] = i;
//				empty_count--;
//				lock.unlock();
//
//				if (empty_count == 0) {
//					return true; // Puzzle solved, terminate recursion
//				}
//
//				if (solve_sudoku(empty_count))
//					return true;
//
//				lock.lock();
//				board[row][col] = 0;
//				empty_count++;
//				lock.unlock();
//
//			}
//		}
//		return false;
//	}
//
//public:
//	SudokuSolverParallel(const string& filename) : SudokuSolver(filename) {}
//bool solve() override {
//
//		//int num_threads = 3;
//		//if (emptyCells > 50) {
//			int num_threads = min(emptyCells / 10+1, static_cast<int>(thread::hardware_concurrency()));
//		//}
//
//		//int num_threads = min(empty_count/10+1, static_cast<int>(thread::hardware_concurrency()));
//		vector<future<bool>> futures;
//
//		for (int i = 0; i < num_threads; i++) {
//			futures.push_back(async(launch::async, [this]() {
//				return solve_sudoku(emptyCells);
//				}));
//		}
//
//		for (auto& f : futures) {
//			if (f.get())
//				return true;
//		}
//		return false;
//	}
//};

// Derived class for solving Sudoku using parallel approach
class SudokuSolverParallel : public SudokuSolver {
private:
	vector<Position> empty_cells;  // List of empty cells on the board
	vector<vector<bool>> row_used;  // Array to track numbers used in each row
	vector<vector<bool>> col_used;  // Array to track numbers used in each column
	vector<vector<vector<bool>>> box_used;  // Array to track numbers used in each box
	mutex board_mutex;  // Mutex for synchronizing access to the board

	bool is_valid(int row, int col, int num) {
		// Check if the number is valid in the given position
		return !row_used[row][num] && !col_used[col][num] && !box_used[row / boxSize][col / boxSize][num];
	}

	bool solve_sudoku(int& empty_count) {
		// Solve the Sudoku using parallel approach
		if (empty_count == 0)
			return true;

		int row = empty_cells[empty_count - 1].first;
		int col = empty_cells[empty_count - 1].second;

		for (int i = 1; i <= boardSize; i++) {
			{
				unique_lock<mutex> lock(board_mutex);
				if (!is_valid(row, col, i))
					continue;
			}

			{
				unique_lock<mutex> lock(board_mutex);
				board[row][col] = i;
				row_used[row][i] = true;
				col_used[col][i] = true;
				box_used[row / boxSize][col / boxSize][i] = true;
				empty_count--;
			}

			if (solve_sudoku(empty_count))
				return true;

			{
				unique_lock<mutex> lock(board_mutex);
				board[row][col] = 0;
				row_used[row][i] = false;
				col_used[col][i] = false;
				box_used[row / boxSize][col / boxSize][i] = false;
				empty_count++;
			}
		}

		return false;
	}

public:
	SudokuSolverParallel(const string& filename) : SudokuSolver(filename) {
		// Initialize data structures for parallel solving
		row_used.resize(boardSize, vector<bool>(boardSize + 1, false));
		col_used.resize(boardSize, vector<bool>(boardSize + 1, false));
		box_used.resize(boxSize, vector<vector<bool>>(boxSize, vector<bool>(boardSize + 1, false)));

		for (int i = 0; i < boardSize; i++) {
			for (int j = 0; j < boardSize; j++) {
				if (board[i][j] == 0) {
					empty_cells.emplace_back(i, j);
				}
				else {
					row_used[i][board[i][j]] = true;
					col_used[j][board[i][j]] = true;
					box_used[i / boxSize][j / boxSize][board[i][j]] = true;
				}
			}
		}
	}

	bool solve() override {
		int num_threads = min(emptyCells / 10 + 1, static_cast<int>(thread::hardware_concurrency()));
		vector<future<bool>> futures;

		for (int i = 0; i < num_threads; i++) {
			futures.push_back(async(launch::async, [this]() {
				int empty_count = emptyCells;
				return solve_sudoku(empty_count);
				}));
		}

		for (auto& f : futures) {
			if (f.get())
				return true;
		}

		return false;
	}
};

void runSolver(SudokuSolver& solver, const string& algorithmName, const string& outputFile) {
	cout << "\n" << "****************************** INPUT GRID - " << algorithmName << " *****************************" << "\n\n";
	solver.printSolution();
	cout << "\n" << "**************************************************************************************" << "\n";

	auto start = chrono::high_resolution_clock::now();
	bool solved = solver.solve();
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

	if (solved) {
		cout << "\n" << "SOLVED!" << "\n";
		cout << "Time taken: " << duration.count() << " microseconds" << endl;
		cout << "****************************** OUTPUT GRID - " << algorithmName << " ****************************" << "\n\n";
		solver.printSolution();
		solver.writeOutput(outputFile, duration.count());
		cout << "\n" << "**************************************************************************************" << "\n";
	}
	else {
		cout << "\n" << "NO SOLUTION FOUND!" << "\n";
	}
}

int main(int argc, char** argv) {
	cout << "\n" << "Sudoku Solver" << "\n" << "developed by Anagha Dhekne" << "\n\n\n";

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <PATH_TO_INPUT_FILE>" << "\n";
		cerr << "Please try again." << "\n";
		exit(-1);
	}

	if (argc != 3) {
		cerr << "Usage: " << argv[0] << " <PATH_TO_INPUT_FILE> <MODE>" << "\n";
		cerr << "        1. <MODE>: " << "\n";
		cerr << "            - 1: sequential mode with backtracking algorithm" << "\n";
		cerr << "            - 2: sequential mode with brute force algorithm" << "\n";
		cerr << "            - 3: parallel mode" << "\n";
		cerr << "            - 4: All" << "\n";
		cerr << "Please try again." << "\n";
		exit(-1);
	}

	int mode = stoi(argv[2]);
	string inputFile = argv[1];
	string outputDirectory = "output/";

	switch (mode) {
	case 1: {
		string outputFile = outputDirectory + "solution_backtracking_" + inputFile.substr(inputFile.rfind('/') + 1);
		SudokuSolverBackTracking solver(inputFile);
		runSolver(solver, "BACKTRACKING", outputFile);
		break;
	}
	case 2: {
		string outputFile = outputDirectory + "solution_bruteforce_" + inputFile.substr(inputFile.rfind('/') + 1);;
		SudokuSolverBruteForce solver(inputFile);
		runSolver(solver, "BRUTEFORCE", outputFile);
		break;
	}
	case 3: {
		string outputFile = outputDirectory + "solution_parallel_" + inputFile.substr(inputFile.rfind('/') + 1);;
		// Parallel mode implementation
		SudokuSolverParallel solver(inputFile);
		runSolver(solver, "PARALLEL", outputFile);
		break;
	}
	case 4: {
		string outputFileBacktracking = outputDirectory + "solution_backtracking_" + inputFile.substr(inputFile.rfind('/') + 1);;
		string outputFileBruteForce = outputDirectory + "solution_bruteforce_" + inputFile.substr(inputFile.rfind('/') + 1);;
		string outputFileParallel = outputDirectory + "solution_parallel_" + inputFile.substr(inputFile.rfind('/') + 1);;

		SudokuSolverBackTracking solverBacktracking(inputFile);
		SudokuSolverBruteForce solverBruteForce(inputFile);
		SudokuSolverParallel solver(inputFile);

		runSolver(solverBacktracking, "BACKTRACKING", outputFileBacktracking);
		runSolver(solverBruteForce, "BRUTEFORCE", outputFileBruteForce);
		runSolver(solver, "PARALLEL", outputFileParallel);

		break;
	}
	default: {
		cerr << "Available options for <MODE>: " << "\n";
		cerr << "            - 1: sequential mode with backtracking algorithm" << "\n";
		cerr << "            - 2: sequential mode with brute force algorithm" << "\n";
		cerr << "            - 3: parallel mode" << "\n";
		cerr << "            - 4: All" << "\n";
		cerr << "Please try again." << "\n";
		exit(-1);
	}
	}

	return 0;
}
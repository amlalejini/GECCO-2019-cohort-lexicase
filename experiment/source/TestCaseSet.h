/**
*   @note The original version of this class was written by Emily Dolson and can be found here: https://github.com/emilydolson/map-elites-gp/blob/master/source/TestcaseSet.h
*/

#ifndef TEST_CASE_SET_H
#define TEST_CASE_SET_H

#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>

#include "base/array.h"
#include "base/vector.h"
#include "tools/string_utils.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

#include "parser.hpp"

template <typename INPUT_TYPE, typename OUTPUT_TYPE>
class TestCaseSet {
protected:
    using input_t = INPUT_TYPE;
    using output_t = OUTPUT_TYPE;
    using test_case_t = std::pair<input_t, output_t>;

    using fun_load_test_case_from_line_str_t = std::function<test_case_t(const std::string &)>;
    using fun_load_test_case_from_line_vec_t = std::function<test_case_t(const emp::vector<std::string> &)>;
    
    emp::vector<test_case_t> test_cases;

    fun_load_test_case_from_line_str_t fun_load_test_case_from_line_str;
    fun_load_test_case_from_line_vec_t fun_load_test_case_from_line_vec;

public:
    // TestCaseSet(const load_test_case_fun_t & load_fun, const std::string & filename) {
    //     fun_load_test_case = load_fun;
    //     LoadTestCases(filename);
    // }

    TestCaseSet(const fun_load_test_case_from_line_str_t & load_fun) {
        fun_load_test_case_from_line_str = load_fun;
        fun_load_test_case_from_line_vec = [](const emp::vector<std::string> &) -> std::pair<input_t, output_t> { return {input_t(), output_t()}; };
    }

    TestCaseSet(const fun_load_test_case_from_line_vec_t & load_fun) {
        fun_load_test_case_from_line_str = [](const std::string &) -> std::pair<input_t, output_t> { return {input_t(), output_t()}; };
        fun_load_test_case_from_line_vec = load_fun;
    }

    TestCaseSet() {
        fun_load_test_case_from_line_str = [](const std::string &) -> std::pair<input_t, output_t> { return {input_t(), output_t()}; };
        fun_load_test_case_from_line_vec = [](const emp::vector<std::string> &) -> std::pair<input_t, output_t> { return {input_t(), output_t()}; };
    }

    test_case_t & operator[](size_t id) { return test_cases[id]; }
    const test_case_t & operator[](size_t id) const { return test_cases[id]; }
    
    size_t GetSize() const { return test_cases.size(); }
    
    /// Get input for testcase given by testID.
    input_t & GetInput(size_t testID) {
        emp_assert(testID < test_cases.size());
        return test_cases[testID].first;
    }
    /// Get input for testcase given by testID.
    const input_t & GetInput(size_t testID) const {
        emp_assert(testID < test_cases.size());
        return test_cases[testID].first;
    }

    /// Get output for test case given by testID.
    output_t & GetOutput(size_t testID) {
        emp_assert(testID < test_cases.size());
        return test_cases[testID].second;
    }
    /// Get output for test case given by testID.
    const output_t & GetOutput(size_t testID) const {
        emp_assert(testID < test_cases.size());
        return test_cases[testID].second;
    }

    /// Get test case set
    emp::vector<test_case_t> & GetTestCaseSet() { return test_cases; }

    void SetLoadFun(const fun_load_test_case_from_line_str_t & load_fun) { fun_load_test_case_from_line_str = load_fun; }
    void SetLoadFun(const fun_load_test_case_from_line_vec_t & load_fun) { fun_load_test_case_from_line_vec = load_fun; }
    
    /// NOTE - in future, deprecate this way of reading things in.
    void LoadTestCases(std::string filename) {
        std::ifstream infile(filename);
        std::string line;
        if (!infile.is_open()) {
            std::cout << "ERROR: " << filename << " did not open correctly." << std::endl;
            return;
        }
        // Ignore header
        getline(infile, line);
        while (getline(infile, line)) {
            test_cases.emplace_back(fun_load_test_case_from_line_str(line));
        }
        infile.close();
    }

    /// NOTE - in future, move forward with this way of reading test case input!
    void LoadTestCasesWithCSVReader(std::string filename) {
        std::ifstream infile(filename);
        aria::csv::CsvParser parser(infile);
        
        bool header = true;
        for (auto & row : parser) {
            if (header) { header = false; continue; } // Skip over header row
            emp::vector<std::string> fields;
            for (auto & field : row) {
                fields.emplace_back(field);
            }
            test_cases.emplace_back(fun_load_test_case_from_line_vec(fields));
        }
    }

    bool EvaluateOnTest(size_t testID, const output_t & out) {
        emp_assert(testID < test_cases.size());
        return out == GetOutput(testID);
    }

};

#endif

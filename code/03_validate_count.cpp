#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <algorithm>
#include <regex>
#include <numeric>
#include <type_traits>

using namespace std;

// stringを指定したdelimiterで分割
vector<string> split(const string &s, char delim){
	vector<string> elems;
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) elems.push_back(item);
	return elems;
}

// stringを指定したdelimiterで分割し、intに変換
vector<int> split_to_int(const string &s, char delim){
	vector<string> parts = split(s, delim);
	vector<int> result;
	bool error_flag = false;
	for(const string &p : parts){
		try{
			result.push_back(stoi(p));
		}catch(const invalid_argument &e){
			error_flag = true;
			// cerr << "Warning: invalid argument in stoi(): \"" << p << "\" in \"" << s << "\"" << endl;
			result.push_back(-1);
		}catch(const out_of_range &e){
			error_flag = true;
			// cerr << "Warning: out of range in stoi(): \"" << p << "\" in \"" << s << "\"" << endl;
			result.push_back(-1);
		}
	}
	if(error_flag){
		// cerr << "Warning [split_to_int]: some parts could not be converted to int in \"" << s << "\"" << endl;
	}

	return result;
}

// TRGT_file クラス
class TRGT_file {
public:
	vector<string> trid;
	vector<vector<string>> motifs;
	vector<vector<int>> coverage;
	vector<vector<string>> mc;
	set<int> AL_set;

	TRGT_file(const string &vcf_path){
		ifstream fin(vcf_path);
		string line;
		while(getline(fin, line)){
			if(line[0] == '#') continue;
			// cout << "Processing line: " << line.substr(0, 20) << endl;
			auto fields = split(line, '\t');

			// [Debugging output] 各フィールドの内容を確認
			// for(size_t i = 0; i < fields.size(); i++){
			// 	cout << "Field " << i << ": " << fields[i] << endl;
			// }
			
			string chrom = fields[0];
			string pos = fields[1];
			string info = fields[7];
			string format = fields[8];
			string sample = fields[9];

			auto info_parts = split(info, ';');
			string trid_field = split(info_parts[0], '=')[1];
			vector<string> motifs_units = split(split(info_parts[2], '=')[1], ',');
			for(const string &motif : motifs_units) AL_set.insert(motif.length());

			auto sample_fields = split(sample, ':');
			auto covs = split_to_int(sample_fields[3], ',');
			auto mcs = split(sample_fields[4], ',');

			trid.push_back(trid_field);
			motifs.push_back(motifs_units);
			coverage.push_back(covs);
			mc.push_back(mcs);
		}
	}
};

// TRGTdenovo_file クラス
class TRGTdenovo_file {
public:
	vector<string> trid;
	vector<int> gt_index, denovo_coverage, father_overlap_coverage, mother_overlap_coverage;
	vector<double> child_ratio;
	vector<string> allele_origin, denovo_status;
	vector<vector<int>> motifs;

	vector<vector<int>> coverage_father, coverage_mother, coverage_child;
	vector<vector<int>> AL_father, AL_mother, AL_child;
	vector<vector<string>> MC_father, MC_mother, MC_child;

	vector<vector<bool>> basic_filtering;
	vector<bool> is_denovo;

	vector<int> unit_count_rows;
	map<int, int> unit_count_rows_map;
	vector<vector<int>> unit_count, unit_count_dup, unit_count_cont;

	TRGTdenovo_file(const string &denovo_path){
		ifstream fin(denovo_path);
		string line;
		getline(fin, line); // skip header
		while(getline(fin, line)){
			auto fields = split(line, '\t');

			// [Debugging output] 各フィールドの内容を確認
			// for(auto &field : fields){
			// 	cout << "Field: " << field << endl;
			// }

			trid.push_back(fields[4]);	// trid
			gt_index.push_back(stoi(fields[23]));	// index
			denovo_coverage.push_back(stoi(fields[6]));
			child_ratio.push_back(stod(fields[10]));
			father_overlap_coverage.push_back(stoi(fields[13]));
			mother_overlap_coverage.push_back(stoi(fields[14]));
			allele_origin.push_back(fields[15]);
			denovo_status.push_back(fields[16]);

			AL_father.push_back(split_to_int(fields[27], ','));
			AL_mother.push_back(split_to_int(fields[28], ','));
			AL_child.push_back(split_to_int(fields[29], ','));

			MC_father.push_back(split(fields[24], ','));
			MC_mother.push_back(split(fields[25], ','));
			MC_child.push_back(split(fields[26], ','));

			basic_filtering.push_back(vector<bool>(5, false));
			is_denovo.push_back(false);
		}
	}

	void get_coverage(const TRGT_file &father, const TRGT_file &mother, const TRGT_file &child, const map<string, int> &trid_map){
		for(const auto &trid_val : trid){
			int idx = trid_map.at(trid_val);
			coverage_father.push_back(father.coverage[idx]);
			coverage_mother.push_back(mother.coverage[idx]);
			coverage_child.push_back(child.coverage[idx]);
		}
	}

	void get_motifs(const TRGT_file &child, const map<string, int> &trid_map){
		for(const auto &trid_val : trid){
			int idx = trid_map.at(trid_val);
			vector<int> unit_lens;
			for(const auto &unit : child.motifs[idx]) unit_lens.push_back(unit.length());
			motifs.push_back(unit_lens);
		}
	}

	void add_statistics_line(){
		// trid, string型
		trid.push_back("Statistics");
		allele_origin.push_back(".");
		denovo_status.push_back(".");

		// double, int, vector<int> 型
		child_ratio.push_back(-1.0);
		gt_index.push_back(-1);
		motifs.push_back({-1});
		for(auto &v : {&denovo_coverage, &father_overlap_coverage, &mother_overlap_coverage}){
			v->push_back(-1);
		}
		for(auto &v : {&coverage_father, &coverage_mother, &coverage_child, &AL_father, &AL_mother, &AL_child}){
			v->push_back({-1, -1});
		}

		// vector<string> 型
		for(auto &v : {&MC_father, &MC_mother, &MC_child}){
			v->push_back({".", "."});
		}

		// bool, vector<bool> 型
		is_denovo.push_back(false);
		basic_filtering.push_back(vector<bool>(5, false));

		// count関係行はcount_unit_size()内で処理する
	}

	void count_unit_size(const TRGT_file &father, const TRGT_file &mother, const TRGT_file &child){
		// unit長の一覧を取得
		set<int> all_units = father.AL_set;
		all_units.insert(mother.AL_set.begin(), mother.AL_set.end());
		all_units.insert(child.AL_set.begin(), child.AL_set.end());
		unit_count_rows.assign(all_units.begin(), all_units.end());
		for(size_t i = 0; i < unit_count_rows.size(); i++) unit_count_rows_map[unit_count_rows[i]] = i;
		for(auto &v : {&unit_count, &unit_count_dup, &unit_count_cont}){
			v->resize(trid.size(), vector<int>(unit_count_rows.size(), 0));
		}

		// [Debugging output] unit_count_rowsの内容を確認
		// for(int i : unit_count_rows){
		// 	cout << i << " ";
		// }cout << endl;

		// [Debugging output] mapの内容を確認
		// cout << "unit_count_rows_map: ";
		// for(const auto &pair : unit_count_rows_map){
		// 	cout << pair.first << "->" << pair.second << " ";
		// }
		// cout << endl;

		for(size_t i = 0; i < trid.size()-1; i++){
			// allele_originが不明の場合，denovo_statusが定まらないのでskip
			if(allele_origin[i] == "" || allele_origin[i] == ".") continue;
			if(allele_origin[i].find("?") != string::npos) continue;		// "?", "F:?", "M:?" の場合対象外
			// filter[0]がfalseの場合も対象外 (2025/07/20追加)
			if(!basic_filtering[i][0]) continue;
			// 例外処理
			if(trid[i] == "Statistics") continue;

			// 対象アレル由来をに対応する親アレルを特定
			auto origin_str = (allele_origin[i].substr(0, 1) == "F") ? MC_father[i][stoi(allele_origin[i].substr(2)) - 1] : MC_mother[i][stoi(allele_origin[i].substr(2)) - 1];
			auto origin_vec = split_to_int(origin_str, '_');

			// origin_vecの検証（例外処理）
			if((int)origin_vec.size() == 0 || origin_vec[0] < 0){
				cerr << "[Warning] (count_unit_size): origin vec is empty for trid " << trid[i] << ", Origin: " << origin_str << endl;
				continue;
			}
			if(origin_vec.size() != motifs[i].size()){
				cerr << "[Error] (count_unit_size): Mismatched origin_vec and motifs size for trid " << trid[i] << ", Origin: " << origin_str << ", Motifs: ";
				for(int m : motifs[i]) cerr << m << " ";
				cerr << endl;
				continue;
			}

			// denovo_vecの検証（例外処理）
			auto denovo_vec = split_to_int(MC_child[i][gt_index[i]], '_');
			if((int)denovo_vec.size() == 0 || denovo_vec[0] < 0){
				cerr << "[Warning] (count_unit_size): denovo vec is empty for trid " << trid[i] << ", Origin: " << origin_str << ", Denovo: " << MC_child[i][gt_index[i]] << endl;
				continue;
			}
			if(origin_vec.size() != denovo_vec.size()){
				cerr << "[Error] (count_unit_size): Mismatched MC lengths for trid " << trid[i] << ", Origin: " << origin_str << ", Denovo: " << MC_child[i][gt_index[i]] << endl;
				continue;
			}

			// 対象アレル数（母数）をカウント
			gt_index.back()++;

			const auto &unit_lens = motifs[i];
			for(size_t j = 0; j < unit_lens.size(); j++){
				int idx = unit_count_rows_map.at(unit_lens[j]);
				unit_count[i][idx] += origin_vec[j];
				unit_count.back()[idx] += origin_vec[j];
			}

			// [Debugging output] 各tridの処理内容を確認
			// cout << "trid; " << trid[i] << ", origin_vec: ";
			// for(int j : origin_vec) cout << j << " ";
			// cout << ", motifs: ";
			// for(int j : motifs[i]) cout << j << " ";
			// cout << endl;

			// 以下de novoアレルの処理
			if(!is_denovo[i]) continue;
			
			for(size_t j = 0; j < unit_lens.size(); j++){
				// 例外処理
				if(denovo_vec[j] < 0){
					cerr << "[Error] (count_unit_size): Negative values found in denovo_vec for trid " << trid[i] << ", denovo_str: \"" << MC_child[i][gt_index[i]] << "\", denovo_vec: ";
					for(int v : denovo_vec) cerr << v << " ";
					cerr << endl;
					continue;
				}
				if(origin_vec[j] < 0){
					cerr << "[Error] (count_unit_size): Negative values found in origin_vec for trid " << trid[i] <<", origin_str: \"" << origin_str << "\", origin_vec: ";
					for(int v : origin_vec) cerr << v << " ";
					cerr << endl;
					continue;
				}

				int diff = denovo_vec[j] - origin_vec[j];
				int idx = unit_count_rows_map.at(unit_lens[j]);
				if(idx < 0 || idx >= (int)unit_count_rows.size()){
					cerr << "[Error] (count_unit_size): Invalid index " << idx << " for unit length " << unit_lens[j] << " in trid " << trid[i] << endl;
					continue;
				}

				// 変異の記録
				if(diff > 0){
					unit_count_dup[i][idx] += diff;
					unit_count_dup.back()[idx] += diff;
					cout << "[INFO] (duplication): trid " << trid[i] << ", unit idx: " << j << ", unit length " << unit_lens[j] << ", count increased by " << diff << endl;
				}
				if(diff < 0){
					unit_count_cont[i][idx] -= diff;
					unit_count_cont.back()[idx] -= diff;
					cout << "[INFO] (contraction): trid " << trid[i] << ", unit idx: " << j << ", unit length " << unit_lens[j] << ", count decreased by " << -diff << endl;
				}
			}
		}
	}
};

void create_trid_map(map<string, int> &trid_map, const TRGT_file &trgt){
	for(size_t i = 0; i < trgt.trid.size(); i++) trid_map[trgt.trid[i]] = i;
}

void filter_denovo(TRGTdenovo_file &d, const map<string, int> &trid_map){
	for(size_t i = 0; i < d.trid.size(); i++){
		auto sum = [](const vector<int> &v){ return accumulate(v.begin(), v.end(), 0); };
		d.basic_filtering[i][0] = sum(d.coverage_father[i]) >= 10 && sum(d.coverage_mother[i]) >= 10 && sum(d.coverage_child[i]) >= 10;

		int c = d.gt_index[i];
		d.basic_filtering[i][1] = !(
			   (find(d.AL_father[i].begin(), d.AL_father[i].end(), d.AL_child[i][c]) != d.AL_father[i].end()
				&& find(d.AL_mother[i].begin(), d.AL_mother[i].end(), d.AL_child[i][1 - c]) != d.AL_mother[i].end())
			|| (find(d.AL_mother[i].begin(), d.AL_mother[i].end(), d.AL_child[i][c]) != d.AL_mother[i].end()
				&& find(d.AL_father[i].begin(), d.AL_father[i].end(), d.AL_child[i][1 - c]) != d.AL_father[i].end())
		);

		d.basic_filtering[i][2] = (d.denovo_status[i] == "Y:+" || d.denovo_status[i] == "Y:-");
		d.basic_filtering[i][3] = (d.denovo_coverage[i] >= 2 && d.child_ratio[i] >= 0.2);

		double father_ratio = (double)d.father_overlap_coverage[i] / max(sum(d.coverage_father[i]), 1);
		double mother_ratio = (double)d.mother_overlap_coverage[i] / max(sum(d.coverage_mother[i]), 1);
		d.basic_filtering[i][4] = (father_ratio < 0.05 && mother_ratio < 0.05);

		d.is_denovo[i] = all_of(d.basic_filtering[i].begin(), d.basic_filtering[i].end(), [](bool x){ return x; });
	}
}

void write_trgtdenovo_to_tsv(const TRGTdenovo_file &d, const string &out_path){
	ofstream fout(out_path);
	size_t N = d.trid.size();
	// 列名の出力
	fout << "trid\tgt_index\tdenovo_coverage\tchild_ratio\tfather_overlap_coverage\tmother_overlap_coverage"
		 << "\tallele_origin\tdenovo_status"
		 << "\tmotifs"
		 << "\tcoverage_father.1\tcoverage_father.2"
		 << "\tcoverage_mother.1\tcoverage_mother.2"
		 << "\tcoverage_child.1\tcoverage_child.2"
		 << "\tAL_father.1\tAL_father.2"
		 << "\tAL_mother.1\tAL_mother.2"
		 << "\tAL_child.1\tAL_child.2"
		 << "\tMC_father.1\tMC_father.2"
		 << "\tMC_mother.1\tMC_mother.2"
		 << "\tMC_child.1\tMC_child.2"
		 << "\tbasic_filtering.1\tbasic_filtering.2\tbasic_filtering.3\tbasic_filtering.4\tbasic_filtering.5"
		 << "\tis_denovo";

	// unit_count 系は unit 長を列名にする
	for(int len : d.unit_count_rows) fout << "\t" << len;
	for(int len : d.unit_count_rows) fout << "\tcont." << len;
	for(int len : d.unit_count_rows) fout << "\tdup." << len;
	fout << "" << endl;

	// 各行を出力
	for(size_t i = 0; i < N; i++){
		// cout << "Outputting row " << i + 1 << " of " << N << ": trid=" << d.trid[i] << endl;
		auto write_int = [&](int value){
			if(value >= 0){fout << "\t" << value;}
			else{
				// cerr << "Warning [output]: negative value for trid " << d.trid[i] << ", value=" << value << endl;
				fout << "\t" << ".";
			}
		};
		auto write_double = [&](double value){
			if(value >= 0){fout << "\t" << value;}
			else{
				// cerr << "Warning [output]: negative value for trid " << d.trid[i] << ", value=" << value << endl;
				fout << "\t" << ".";
			}
		};

		fout << d.trid[i];
		write_int(d.gt_index[i]);
		write_int(d.denovo_coverage[i]);
		write_double(d.child_ratio[i]);
		write_int(d.father_overlap_coverage[i]);
		write_int(d.mother_overlap_coverage[i]);
		fout << "\t" << d.allele_origin[i] << "\t" << d.denovo_status[i];

		auto write_vec = [&](const auto &vec){
			fout << "\t";
			for(int j = 0; j < (int)vec[i].size(); j++){
				fout << vec[i][j];
				if(j != (int)vec[i].size() - 1) fout << ",";
			}
		};
		auto write_pair = [&](const auto &vec, const string &name){
			if((int)vec[i].size() == 0){
				// cerr << "Warning [output]: vector size is zero for trid " << d.trid[i] << ", field=" << name << endl;
				fout << "\t.\t.";
				return;
			}
			using T = decay_t<decltype(vec[i][0])>;
			if((int)vec[i].size() == 1){
				// cerr << "Warning [output]: vector size is 1 for trid " << d.trid[i] << ", field=" << name << ", vec[0]=" << vec[i][0] << endl;
				if constexpr (is_integral_v<T>){
					if(vec[i][0] >= 0){
						fout << "\t" << vec[i][0] << "\t.";
					}else{
						fout << "\t.\t.";
					}
				}else{
					fout << "\t" << vec[i][0] << "\t.";
				}
				return;
			}
			if constexpr (is_integral_v<T>){
				if(vec[i][0] >= 0){
					fout << "\t" << vec[i][0];
				}else{
					fout << "\t.";
				}
			}else{
				fout << "\t" << vec[i][0];
			}
			if constexpr (is_integral_v<T>){
				if(vec[i][1] >= 0){
					fout << "\t" << vec[i][1];
				}else{
					fout << "\t.";
				}
			}else{
				fout << "\t" << vec[i][1];
			}
		};

		write_vec(d.motifs);
		write_pair(d.coverage_father, "coverage_father");
		write_pair(d.coverage_mother, "coverage_mother");
		write_pair(d.coverage_child, "coverage_child");
		write_pair(d.AL_father, "AL_father");
		write_pair(d.AL_mother, "AL_mother");
		write_pair(d.AL_child, "AL_child");
		write_pair(d.MC_father, "MC_father");
		write_pair(d.MC_mother, "MC_mother");
		write_pair(d.MC_child, "MC_child");

		for(int k = 0; k < 5; k++) fout << "\t" << d.basic_filtering[i][k];
		fout << "\t" << d.is_denovo[i];

		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count[i][j];
		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count_cont[i][j];
		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count_dup[i][j];
		fout << endl;
	}

	fout.close();
}

void write_tsv_count_only(const TRGTdenovo_file &d, const string &out_path){
	ofstream fout(out_path);
	size_t N = d.trid.size();
	// 列名の出力
	fout << "trid\t";

	// unit_count 系は unit 長を列名にする
	for(int len : d.unit_count_rows) fout << "\t" << len;
	for(int len : d.unit_count_rows) fout << "\tcont." << len;
	for(int len : d.unit_count_rows) fout << "\tdup." << len;

	fout << "\tbasic_filtering.1\tbasic_filtering.2\tbasic_filtering.3\tbasic_filtering.4\tbasic_filtering.5"
		 << "\tis_denovo";

	fout << "" << endl;

	// 各行を出力
	for(size_t i = 0; i < N; i++){
		fout << d.trid[i] << "\t";

		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count[i][j];
		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count_cont[i][j];
		for(int j = 0; j < (int)d.unit_count_rows.size(); j++) fout << "\t" << d.unit_count_dup[i][j];
		
		for(int k = 0; k < 5; k++) fout << "\t" << d.basic_filtering[i][k];
		fout << "\t" << d.is_denovo[i];
		
		fout << endl;
	}

	fout.close();
}

int main(int argc, char *argv[]){
	if(argc != 5){
		cerr << "Usage: " << argv[0] << " <trgt_denovo_file> <child_file> <father_file> <mother_file>" << endl;
		return 1;
	}

	cout << "Starting TRGTdenovo_file processing..." << endl;

	string trgt_denovo_file = argv[1], child_file = argv[2], father_file = argv[3], mother_file = argv[4];
	TRGT_file child(child_file), father(father_file), mother(mother_file);
	cout << "Created TRGT_file objects" << endl;

	map<string, int> trid_map;
	create_trid_map(trid_map, child);

	TRGTdenovo_file denovo(trgt_denovo_file);
	cout << "Created TRGTdenovo_file object" << endl;
	denovo.get_coverage(father, mother, child, trid_map);
	denovo.get_motifs(child, trid_map);
	cout << "Got coverage and structure" << endl;
	filter_denovo(denovo, trid_map);
	cout << "Filtering completed" << endl;
	denovo.add_statistics_line();
	cout << "Added statistics line" << endl;

	denovo.count_unit_size(father, mother, child);
	cout << "Counted all unit sizes" << endl;

	write_trgtdenovo_to_tsv(denovo, "validate_count.tsv");
	write_tsv_count_only(denovo, "result_only.tsv");
}

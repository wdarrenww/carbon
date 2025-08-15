#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../core/system.hpp"
#include "../core/elements.hpp"

namespace carbon {

class XYZReader {
public:
    static bool read_file(const std::string& filename, System& system) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        system.atoms.clear();
        system.bonds.clear();
        
        std::string line;
        int atom_count = 0;
        
        // read number of atoms
        if (!std::getline(file, line)) {
            return false;
        }
        
        std::stringstream ss(line);
        if (!(ss >> atom_count) || atom_count <= 0) {
            return false;
        }
        
        // skip comment line
        if (!std::getline(file, line)) {
            return false;
        }
        
        // read atoms
        const ElementTable& table = get_element_table();
        
        for (int i = 0; i < atom_count; ++i) {
            if (!std::getline(file, line)) {
                return false;
            }
            
            std::stringstream atom_ss(line);
            std::string element;
            double x, y, z;
            
            if (!(atom_ss >> element >> x >> y >> z)) {
                return false;
            }
            
            if (!table.has_element(element)) {
                continue;  // skip unknown elements
            }
            
            Atom atom = table.create_atom(element, Vec3(x, y, z));
            system.add_atom(atom);
        }
        
        file.close();
        return true;
    }
    
    static bool read_frame(std::ifstream& file, System& system, bool append = false) {
        if (!append) {
            system.atoms.clear();
        }
        
        std::string line;
        int atom_count = 0;
        
        // read number of atoms
        if (!std::getline(file, line)) {
            return false;
        }
        
        std::stringstream ss(line);
        if (!(ss >> atom_count) || atom_count <= 0) {
            return false;
        }
        
        // skip comment line
        if (!std::getline(file, line)) {
            return false;
        }
        
        // read atoms
        const ElementTable& table = get_element_table();
        
        for (int i = 0; i < atom_count; ++i) {
            if (!std::getline(file, line)) {
                return false;
            }
            
            std::stringstream atom_ss(line);
            std::string element;
            double x, y, z;
            
            if (!(atom_ss >> element >> x >> y >> z)) {
                return false;
            }
            
            if (!table.has_element(element)) {
                continue;
            }
            
            if (append) {
                Atom atom = table.create_atom(element, Vec3(x, y, z));
                system.add_atom(atom);
            } else {
                if (i < system.num_atoms()) {
                    system.atoms[i].r = Vec3(x, y, z);
                }
            }
        }
        
        return true;
    }
};

class XYZWriter {
public:
    static bool write_file(const std::string& filename, const System& system, 
                          const std::string& comment = "") {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        write_frame(file, system, comment);
        file.close();
        return true;
    }
    
    static void write_frame(std::ofstream& file, const System& system, 
                           const std::string& comment = "") {
        const ElementTable& table = get_element_table();
        
        // write number of atoms
        file << system.num_atoms() << "\n";
        
        // write comment line
        if (comment.empty()) {
            file << "carbon molecular simulation\n";
        } else {
            file << comment << "\n";
        }
        
        // write atoms
        file << std::fixed << std::setprecision(6);
        for (const auto& atom : system.atoms) {
            std::string symbol = table.get_symbol(atom.Z);
            if (symbol.empty()) {
                symbol = "X";  // unknown element
            }
            
            file << symbol << " " 
                 << atom.r.x << " " 
                 << atom.r.y << " " 
                 << atom.r.z << "\n";
        }
    }
    
    static bool append_frame(const std::string& filename, const System& system,
                            const std::string& comment = "") {
        std::ofstream file(filename, std::ios::app);
        if (!file.is_open()) {
            return false;
        }
        
        write_frame(file, system, comment);
        file.close();
        return true;
    }
};

class TrajectoryWriter {
private:
    std::ofstream file;
    std::string filename;
    int frame_count;
    
public:
    explicit TrajectoryWriter(const std::string& fname) 
        : filename(fname), frame_count(0) {}
    
    ~TrajectoryWriter() {
        close();
    }
    
    bool open() {
        file.open(filename);
        frame_count = 0;
        return file.is_open();
    }
    
    void close() {
        if (file.is_open()) {
            file.close();
        }
    }
    
    bool is_open() const {
        return file.is_open();
    }
    
    void write_frame(const System& system, double time = 0.0, double energy = 0.0) {
        if (!file.is_open()) return;
        
        std::stringstream comment;
        comment << "frame " << frame_count 
                << ", time=" << std::fixed << std::setprecision(3) << time << " fs";
        if (energy != 0.0) {
            comment << ", energy=" << std::setprecision(6) << energy << " kcal/mol";
        }
        
        XYZWriter::write_frame(file, system, comment.str());
        frame_count++;
        file.flush();
    }
    
    int get_frame_count() const {
        return frame_count;
    }
};

}
#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "../core/system.hpp"
#include "../core/elements.hpp"
#include "xyz_io.hpp"

namespace carbon {

class InteractiveBuilder {
public:
    static void build_system(System& system) {
        std::cout << "Interactive molecular builder\n";
        std::cout << "Commands:\n";
        std::cout << "  add <element> <x> <y> <z>  - add atom at position\n";
        std::cout << "  bond <i> <j> <k> <r0>      - add bond between atoms i and j\n";
        std::cout << "  list                       - list all atoms\n";
        std::cout << "  bonds                      - list all bonds\n";
        std::cout << "  remove <i>                 - remove atom i\n";
        std::cout << "  clear                      - clear all atoms and bonds\n";
        std::cout << "  box <x> <y> <z>           - set box size\n";
        std::cout << "  periodic <on|off>         - toggle periodic boundaries\n";
        std::cout << "  save <filename>           - save to xyz file\n";
        std::cout << "  load <filename>           - load from xyz file\n";
        std::cout << "  done                      - finish building\n";
        std::cout << "  help                      - show this help\n\n";
        
        print_available_elements();
        
        std::string line;
        while (true) {
            std::cout << "> ";
            if (!std::getline(std::cin, line)) {
                break;
            }
            
            if (line.empty()) continue;
            
            std::stringstream ss(line);
            std::string command;
            ss >> command;
            
            if (command == "done" || command == "exit" || command == "quit") {
                break;
            } else if (command == "help") {
                print_help();
            } else if (command == "add") {
                handle_add_atom(ss, system);
            } else if (command == "bond") {
                handle_add_bond(ss, system);
            } else if (command == "list") {
                list_atoms(system);
            } else if (command == "bonds") {
                list_bonds(system);
            } else if (command == "remove") {
                handle_remove_atom(ss, system);
            } else if (command == "clear") {
                system.atoms.clear();
                system.bonds.clear();
                std::cout << "system cleared\n";
            } else if (command == "box") {
                handle_set_box(ss, system);
            } else if (command == "periodic") {
                handle_periodic(ss, system);
            } else if (command == "save") {
                handle_save(ss, system);
            } else if (command == "load") {
                handle_load(ss, system);
            } else {
                std::cout << "unknown command: " << command << "\n";
                std::cout << "type 'help' for available commands\n";
            }
        }
    }

private:
    static void print_help() {
        std::cout << "\nInteractive molecular builder commands:\n";
        std::cout << "  add <element> <x> <y> <z>  - add atom at position (angstroms)\n";
        std::cout << "  bond <i> <j> <k> <r0>      - add bond (i,j=atom indices, k=spring const, r0=eq length)\n";
        std::cout << "  list                       - list all atoms with indices\n";
        std::cout << "  bonds                      - list all bonds\n";
        std::cout << "  remove <i>                 - remove atom at index i\n";
        std::cout << "  clear                      - clear all atoms and bonds\n";
        std::cout << "  box <x> <y> <z>           - set simulation box size (angstroms)\n";
        std::cout << "  periodic <on|off>         - toggle periodic boundary conditions\n";
        std::cout << "  save <filename>           - save current system to xyz file\n";
        std::cout << "  load <filename>           - load system from xyz file\n";
        std::cout << "  done                      - finish building and continue\n";
        std::cout << "  help                      - show this help\n\n";
    }
    
    static void print_available_elements() {
        std::cout << "Available elements: H, C, N, O, Na, Cl\n\n";
    }
    
    static void handle_add_atom(std::stringstream& ss, System& system) {
        std::string element;
        double x, y, z;
        
        if (!(ss >> element >> x >> y >> z)) {
            std::cout << "usage: add <element> <x> <y> <z>\n";
            return;
        }
        
        const ElementTable& table = get_element_table();
        if (!table.has_element(element)) {
            std::cout << "unknown element: " << element << "\n";
            print_available_elements();
            return;
        }
        
        Atom atom = table.create_atom(element, Vec3(x, y, z));
        int index = system.add_atom(atom);
        
        std::cout << "added " << element << " atom at (" << x << ", " << y << ", " << z 
                  << ") with index " << index << "\n";
    }
    
    static void handle_add_bond(std::stringstream& ss, System& system) {
        int i, j;
        double k, r0;
        
        if (!(ss >> i >> j >> k >> r0)) {
            std::cout << "usage: bond <i> <j> <k> <r0>\n";
            std::cout << "  i, j: atom indices\n";
            std::cout << "  k: spring constant (kcal/mol/angstrom^2)\n";
            std::cout << "  r0: equilibrium bond length (angstroms)\n";
            return;
        }
        
        if (i < 0 || i >= system.num_atoms() || j < 0 || j >= system.num_atoms()) {
            std::cout << "invalid atom indices. valid range: 0 to " << (system.num_atoms() - 1) << "\n";
            return;
        }
        
        if (i == j) {
            std::cout << "cannot bond atom to itself\n";
            return;
        }
        
        system.add_bond(i, j, k, r0);
        double current_distance = system.distance_between_atoms(i, j);
        
        std::cout << "added bond between atoms " << i << " and " << j 
                  << " (k=" << k << ", r0=" << r0 << ", current distance=" << current_distance << ")\n";
    }
    
    static void list_atoms(const System& system) {
        if (system.num_atoms() == 0) {
            std::cout << "no atoms in system\n";
            return;
        }
        
        const ElementTable& table = get_element_table();
        std::cout << "atoms in system:\n";
        std::cout << "idx  element  x        y        z        mass\n";
        std::cout << "---  -------  -------  -------  -------  -------\n";
        
        for (int i = 0; i < system.num_atoms(); ++i) {
            const Atom& atom = system.atoms[i];
            std::string symbol = table.get_symbol(atom.Z);
            
            std::cout << std::setw(3) << i << "  "
                      << std::setw(7) << symbol << "  "
                      << std::setw(7) << std::fixed << std::setprecision(3) << atom.r.x << "  "
                      << std::setw(7) << std::fixed << std::setprecision(3) << atom.r.y << "  "
                      << std::setw(7) << std::fixed << std::setprecision(3) << atom.r.z << "  "
                      << std::setw(7) << std::fixed << std::setprecision(3) << atom.mass << "\n";
        }
    }
    
    static void list_bonds(const System& system) {
        if (system.num_bonds() == 0) {
            std::cout << "no bonds in system\n";
            return;
        }
        
        std::cout << "bonds in system:\n";
        std::cout << "i    j    k         r0        current r\n";
        std::cout << "---  ---  --------  --------  ---------\n";
        
        for (const auto& bond : system.bonds) {
            double current_r = system.distance_between_atoms(bond.i, bond.j);
            
            std::cout << std::setw(3) << bond.i << "  "
                      << std::setw(3) << bond.j << "  "
                      << std::setw(8) << std::fixed << std::setprecision(3) << bond.k << "  "
                      << std::setw(8) << std::fixed << std::setprecision(3) << bond.r0 << "  "
                      << std::setw(9) << std::fixed << std::setprecision(3) << current_r << "\n";
        }
    }
    
    static void handle_remove_atom(std::stringstream& ss, System& system) {
        int index;
        if (!(ss >> index)) {
            std::cout << "usage: remove <index>\n";
            return;
        }
        
        if (index < 0 || index >= system.num_atoms()) {
            std::cout << "invalid atom index. valid range: 0 to " << (system.num_atoms() - 1) << "\n";
            return;
        }
        
        system.remove_atom(index);
        std::cout << "removed atom " << index << "\n";
    }
    
    static void handle_set_box(std::stringstream& ss, System& system) {
        double x, y, z;
        if (!(ss >> x >> y >> z)) {
            std::cout << "usage: box <x> <y> <z>\n";
            return;
        }
        
        system.box_size = Vec3(x, y, z);
        std::cout << "box size set to (" << x << ", " << y << ", " << z << ")\n";
    }
    
    static void handle_periodic(std::stringstream& ss, System& system) {
        std::string setting;
        if (!(ss >> setting)) {
            std::cout << "usage: periodic <on|off>\n";
            std::cout << "current setting: " << (system.periodic ? "on" : "off") << "\n";
            return;
        }
        
        if (setting == "on" || setting == "true" || setting == "1") {
            system.periodic = true;
            std::cout << "periodic boundary conditions enabled\n";
        } else if (setting == "off" || setting == "false" || setting == "0") {
            system.periodic = false;
            std::cout << "periodic boundary conditions disabled\n";
        } else {
            std::cout << "usage: periodic <on|off>\n";
        }
    }
    
    static void handle_save(std::stringstream& ss, System& system) {
        std::string filename;
        if (!(ss >> filename)) {
            std::cout << "usage: save <filename>\n";
            return;
        }
        
        if (XYZWriter::write_file(filename, system)) {
            std::cout << "system saved to " << filename << "\n";
        } else {
            std::cout << "failed to save to " << filename << "\n";
        }
    }
    
    static void handle_load(std::stringstream& ss, System& system) {
        std::string filename;
        if (!(ss >> filename)) {
            std::cout << "usage: load <filename>\n";
            return;
        }
        
        if (XYZReader::read_file(filename, system)) {
            std::cout << "system loaded from " << filename << "\n";
            std::cout << "loaded " << system.num_atoms() << " atoms\n";
        } else {
            std::cout << "failed to load from " << filename << "\n";
        }
    }
};

}
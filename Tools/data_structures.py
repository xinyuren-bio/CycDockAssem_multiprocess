import os
import numpy as np

class Fragment:
    def __init__(self, pdb_path):
        self.path = pdb_path
        self.name = os.path.basename(pdb_path)
        self.coords, self.atom_names, self.residue_info, self.key_atom_indices = None, [], {}, {}
        self.atom_count, self.residue_count, self.residues = 0, 0, []
        self._manual_parse_pdb()

    def is_valid(self):
        required_keys = ['N_N', 'CA_N', 'C_N', 'N_C', 'CA_C', 'C_C']
        return (self.coords is not None and self.atom_count > 0 and
                all(key in self.key_atom_indices and self.key_atom_indices[key] is not None for key in required_keys))

    def _manual_parse_pdb(self):
        if not os.path.exists(self.path): return

        coords_list, temp_residues, self.atom_names = [], {}, []
        try:
            with open(self.path, 'r', errors='ignore') as f:
                atom_index_counter = 0
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            atom_name = line[12:16].strip()
                            
                            if atom_name not in ['N', 'CA', 'C', 'O']:
                                continue 

                            coords_list.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                            
                            self.atom_names.append(atom_name)
                            res_name = line[17:20].strip()
                            chain_id = line[21:22].strip()
                            res_seq = int(line[22:26])
                            res_id = (chain_id, res_seq)
                            
                            if res_id not in temp_residues:
                                temp_residues[res_id] = {'res_name': res_name, 'start_index': atom_index_counter, 'atoms': {}}
                            
                            temp_residues[res_id]['atoms'][self.atom_names[-1]] = atom_index_counter
                            temp_residues[res_id]['end_index'] = atom_index_counter
                            atom_index_counter += 1
                        except (ValueError, IndexError): continue
            
            self.atom_count = len(coords_list)
            if self.atom_count == 0: return
            self.coords = np.array(coords_list, dtype=np.float32)
            
            self.residues = sorted(temp_residues.keys(), key=lambda r: r[1])
            self.residue_info = {res_id: info for res_id, info in temp_residues.items()}
            self.residue_count = len(self.residues)
            if not self.residues: return
            n_res_id, c_res_id = self.residues[0], self.residues[-1]
            n_term_atoms, c_term_atoms = self.residue_info[n_res_id]['atoms'], self.residue_info[c_res_id]['atoms']
            self.key_atom_indices = {
                'N_N': n_term_atoms.get('N'), 'CA_N': n_term_atoms.get('CA'), 'C_N': n_term_atoms.get('C'),
                'N_C': c_term_atoms.get('N'), 'CA_C': c_term_atoms.get('CA'), 'C_C': c_term_atoms.get('C')
            }
        except Exception:
            self.coords = None


    def get_core_indices(self, sup_type):
        if sup_type == 0:
            if self.residue_count <= 2: return np.array([], dtype=int)
            start_index = self.residue_info[self.residues[1]]['start_index']
            end_index = self.residue_info[self.residues[-1]]['start_index']
            return np.arange(start_index, end_index)
        else:
            start_idx, end_idx = self.key_atom_indices.get('C_N'), self.key_atom_indices.get('N_C')
            if start_idx is None or end_idx is None: return np.array([], dtype=int)
            return np.arange(start_idx, end_idx + 1)


    def get_p1_clash_indices(self):
        if self.residue_count <= 1: return np.array([], dtype=int)
        start_index = self.residue_info[self.residues[1]]['start_index']
        return np.arange(start_index, self.atom_count)

    def get_p2_clash_indices(self):
        if self.residue_count <= 1: return np.array([], dtype=int)
        end_index = self.residue_info[self.residues[-1]]['start_index']
        return np.arange(0, end_index)

    def get_p1_clash_indices_type2(self):
        if self.residue_count <= 1: return np.array([], dtype=int)
        start_index = self.residue_info[self.residues[1]]['start_index'] + 1
        return np.arange(start_index, self.atom_count)

    def get_p2_clash_indices_type2(self):
        if self.residue_count <= 2: return np.array([], dtype=int)
        
        c_atom_index_to_exclude = self.residue_info[self.residues[-2]]['atoms'].get('C')
        
        end_index = self.residue_info[self.residues[-1]]['start_index']
        all_indices = np.arange(0, end_index)
        
        if c_atom_index_to_exclude is not None:
            return all_indices[all_indices != c_atom_index_to_exclude]
        else:
            return all_indices

    def get_ca_coords(self):
        ca_indices = [i for i, name in enumerate(self.atom_names) if name == 'CA']
        return self.coords[ca_indices] if ca_indices else np.array([])
    
    def get_n_term_sidechain_indices(self):
        if self.residue_count == 0: return np.array([], dtype=int)
        n_term_res_info = self.residue_info[self.residues[0]]
        start_index = n_term_res_info['start_index'] + 4
        end_index = n_term_res_info['end_index']
        return np.arange(start_index, end_index + 1) if start_index <= end_index else np.array([], dtype=int)

    def get_c_term_sidechain_indices(self):
        if self.residue_count == 0: return np.array([], dtype=int)
        c_term_res_info = self.residue_info[self.residues[-1]]
        start_index = c_term_res_info['start_index'] + 4
        return np.arange(start_index, self.atom_count)

class FullAtomFragment:
    def __init__(self, pdb_path):
        self.path = pdb_path
        self.name = os.path.basename(pdb_path)
        self.coords, self.atom_names, self.residue_info = None, [], {}
        self.atom_count, self.residue_count, self.residues = 0, 0, []
        self._manual_parse_pdb()

    def _manual_parse_pdb(self):
        if not os.path.exists(self.path): return

        coords_list, temp_residues, self.atom_names = [], {}, []
        try:
            with open(self.path, 'r', errors='ignore') as f:
                atom_index_counter = 0
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            coords_list.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                            
                            self.atom_names.append(line[12:16].strip())
                            res_name = line[17:20].strip()
                            chain_id = line[21:22].strip()
                            res_seq = int(line[22:26])
                            res_id = (chain_id, res_seq)
                            
                            if res_id not in temp_residues:
                                temp_residues[res_id] = {'res_name': res_name, 'start_index': atom_index_counter, 'atoms': {}}
                            
                            temp_residues[res_id]['atoms'][self.atom_names[-1]] = atom_index_counter
                            temp_residues[res_id]['end_index'] = atom_index_counter
                            atom_index_counter += 1
                        except (ValueError, IndexError): continue
            
            self.atom_count = len(coords_list)
            if self.atom_count > 0:
                self.coords = np.array(coords_list, dtype=np.float32)
                self.residues = sorted(temp_residues.keys(), key=lambda r: r[1])
                self.residue_info = {res_id: info for res_id, info in temp_residues.items()}
                self.residue_count = len(self.residues)
        except Exception:
            self.coords = None

    def get_p1_clash_indices(self):
        if self.residue_count <= 1: return np.array([], dtype=int)
        start_index = self.residue_info[self.residues[1]]['start_index']
        return np.arange(start_index, self.atom_count)

    def get_p2_clash_indices(self):
        if self.residue_count <= 1: return np.array([], dtype=int)
        end_index = self.residue_info[self.residues[-1]]['start_index']
        return np.arange(0, end_index)

    def get_n_term_sidechain_indices(self):
        if self.residue_count == 0: return np.array([], dtype=int)
        n_term_res_info = self.residue_info[self.residues[0]]
        start_index = n_term_res_info['start_index'] + 4
        end_index = n_term_res_info['end_index']
        return np.arange(start_index, end_index + 1) if start_index <= end_index else np.array([], dtype=int)

    def get_c_term_sidechain_indices(self):
        if self.residue_count == 0: return np.array([], dtype=int)
        c_term_res_info = self.residue_info[self.residues[-1]]
        start_index = c_term_res_info['start_index'] + 4
        return np.arange(start_index, self.atom_count)
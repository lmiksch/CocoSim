# -*- coding: utf-8 -*-
"""
This module contains all the different util functions. 
"""
#import RNA

def is_balanced_structure(s):
    stack = []
    for char in s:
        if char == '(':
            stack.append('(')
        elif char == ')':
            if not stack:
                return False  # Unmatched closing parenthesis
            stack.pop()
    
    return len(stack) == 0

def is_legal_dotbracket(structure):
    """
    Check if a dot-bracket structure is legal.
    
    A legal structure only contains the characters '(', ')', and '.',
    and every opening parenthesis '(' must be matched with a closing ')'.

    Args:
        structure (str): The dot-bracket structure to check.

    Returns:
        bool: True if the structure is legal, False otherwise.
    """
    allowed_chars = {'(', ')', '.'}
    
    # Check for any invalid characters
    for char in structure:
        if char not in allowed_chars:
            return False
    
    # Check for balanced parentheses
    stack = []
    for char in structure:
        if char == '(':
            stack.append('(')
        elif char == ')':
            if not stack:
                return False  # Unmatched closing parenthesis
            stack.pop()
    
    # If stack is empty, all opening parentheses were properly closed
    return len(stack) == 0


def split_into_modules(input_string):
    import re
    # Use regular expression to split the string by the spacers
    pattern = re.compile(r'\s*S\d+\s*')
    modules = pattern.split(input_string)
    
    # Remove any leading/trailing whitespace from each module
    modules = [module.strip() for module in modules if module.strip()]
    
    return modules


def structure_to_pairtable(structure):
    """
    Convert a dot-bracket structure into a pair table.

    Args:
        structure (str): Dot-bracket string representation of an RNA secondary structure.

    Returns:
        list: A pair table where:
              - pair_table[0] is the number of bases.
              - For each i >= 1, pair_table[i] is the index of the base paired with i,
                or 0 if the base is unpaired.
                
    Raises:
        ValueError: If the structure contains unbalanced parentheses.
    """
    n = len(structure)
    pair_table = [0] * (n + 1)
    pair_table[0] = n  # Store the length at index 0
    stack = []

    # Iterate over the structure (using 1-indexing for the pair table)
    for i, char in enumerate(structure, start=1):
        if char == '(':
            # Push the position onto the stack
            stack.append(i)
        elif char == ')':
            # Pop the last opening bracket's index and form a pair
            if not stack:
                raise ValueError(f"Invalid structure{structure = }: Unbalanced parentheses (extra closing bracket at position {i = })")
            j = stack.pop()
            pair_table[i] = j
            pair_table[j] = i
        elif char == '.':
            # Unpaired base
            pair_table[i] = 0
        else:
            # If the character is unexpected, raise an error
            raise ValueError("Invalid character '{}' in structure at position {}".format(char, i))
    
    if stack:
        raise ValueError("Invalid structure: Unbalanced parentheses (extra opening bracket(s) remain)")

    return pair_table


def pairtablepath_to_dotbracket(path_pairtable):
    """
    Converts a pairtable path (list of pairtable representations) 
    into a dot-bracket path (list of strings).

    Args:
        path_pairtable (list): List of pairtable representations.
            Each pairtable is a list of integers, where:
                - The first element is the length (n)
                - Elements 1..n indicate the pairing partner (or 0 if unpaired)

    Returns:
        list: A list of dot-bracket strings corresponding to each pairtable.
    """
    dotbracket_path = []
    for pt in path_pairtable:
        # pt[0] holds the number of positions in this structure.
        n = pt[0]
        # Initialize an array to hold the dot-bracket characters.
        structure_chars = [''] * n
        for i in range(1, n+1):
            if pt[i] == 0:
                structure_chars[i-1] = '.'
            elif pt[i] > i:
                structure_chars[i-1] = '('
            else:
                structure_chars[i-1] = ')'
        dotbracket_path.append("".join(structure_chars))
    return dotbracket_path


def path_to_pairtablepath(path):
    """ Function which turns a structural path in dot bracket annotation into a pairtable path

    args: 
        path(list): List of the folding path in dot bracket
    output: 
        pairtable_path(List): List of the folding path using the pairtable annotation

    """
    pairtable_path = [[] for x in path]
    pairtable_path[0] = [1,0]
    for x in range(1,len(path)):
        pairtable_path[x] = [0 for i in range(len(path[x])+1)]
        pairtable_path[x][0] = len(path[x]) 
        for z in range(len(path[x])):
            if path[x] == ".":
                pairtable_path[x][z] = 0
            o = 1
            i = 0
            match = False
            if path[x][z] == "(":
                while match == False:
                    i += 1
                    if path[x][z+i] == "(":
                        o += 1
                    elif path[x][z+i] == ")" and o != 0:
                        o -= 1

                    if path[x][z+i] == ")" and o == 0:
                        
                        pairtable_path[x][z+1] = z + i + 1
                        pairtable_path[x][z+i+1] = z + 1
                        match = True
                    
    return pairtable_path


def find_connected_modules(acfp):
    """ Identifies connected Modules and returns the indices of the connected modules in a list. Each sublist corresponds to connected modules. 


    Args: 
        acfp(list): abstract folding path in the form of a pairtable
    
    Returns: 
        connected_modules(list): sublist corresponds to connected modules with indices 
    """

    all_module = list(range(1,len(acfp)+1))
    connected_modules = []
    curr_module = 1
    i,j = 0,0 
    while all_module:
        all_list = [all_module[0]]
        connected_modules.append(all_list)
        all_module.remove(all_module[0])
        
        while curr_module <= len(acfp):
            for step in range(1,len(acfp)): 
                
                try:
                    if acfp[step][curr_module] != 0 and acfp[step][curr_module] in all_module:
                        connected_modules[i].append(acfp[step][curr_module])
                        all_module.remove(acfp[step][curr_module])
                except:
                    continue
            
            j += 1 

            try:
                curr_module = connected_modules[i][j]
                

            except:
                if len(all_module) == 0:    
                    return connected_modules
                i += 1 
                j = 0 
                all_list = [all_module[0]]
                connected_modules.append(all_list)
                all_module.remove(all_module[0])
                curr_module = connected_modules[i][j]
                
                
    print("Error Something went wrong please check input")


def cv_db2kernel(domain_seq,dot_bracket_seq):
    """
    Takes a dot-bracket sequence and a domain level sequence and converts it into a kernel annotated sequence.

    Args:
        domain_seq(str): domain level seq with space annotation
        dot_bracket_seq(str): dot bracket seq 

    Returns: 
        kernel_seq(str): kernel annotated string with space annotation
    """
    

    
    split_dot_bracket = list(dot_bracket_seq)

    bps = dot_bracket_seq.count('(')
    kernel_annotation = [[] for x in split_dot_bracket]

    open_domains = [] 

    for i, i_domain in enumerate(domain_seq):

        if split_dot_bracket[i] == "(":
            kernel_annotation[i].append(str(i_domain + "("))

        if split_dot_bracket[i] == ")":
            kernel_annotation[i].append(")")
    
        if split_dot_bracket[i] == "." or split_dot_bracket[i] == "x":
            kernel_annotation[i].append(domain_seq[i])
        
    joined_kernel_annotation = []

    


    for x in kernel_annotation:
        joined_kernel_annotation.append(x[0])

    

    return " ".join(joined_kernel_annotation)


def kernel_to_dot_bracket(kernel_seq):
    """
    Takes a kernel annotated sequence and converts it back to a dot-bracket sequence.

    Args:
        kernel_seq(str): kernel annotated string with space annotation

    Returns:
        dot_bracket_seq(str): dot bracket sequence
    """
    kernel_annotations = kernel_seq.split()
    dot_bracket_seq = ""

    for annotation in kernel_annotations:
        if "(" in annotation:
            dot_bracket_seq += "("
        elif ")" in annotation:
            dot_bracket_seq += ")"
        else:
            dot_bracket_seq += "."


    
    return dot_bracket_seq



def is_star_pair(a,b):
    """Checks if two domain pair in Star Annotation"""
    if a[0] == b[0]:
        if len(a) < len(b) or len(a) > len(b):
            return True
        
def only_logic_domain_struct(seq,path):
    """Takes a seq and the extended folding path only returns the folding path of only the logic domains

    Args: 
        seq(list): sequence in form of a list where each entry in the list corresponds to one module
        path(str): module folding path in form of a list where each entry corresponds to one module being transcribed

    """
    logic_domain_struct = []
    for x in range(len(path)):
        
        if seq[x].startswith("L"):
            
            logic_domain_struct.append(path[x])

    
    logic_domain_struct = "".join(logic_domain_struct)
    return logic_domain_struct

def acfp_terminal_input():

    '''Function to let user put in an abstract folding path in the terminal while checking if it's viable
    '''
    acfp = []
    while True:
        print("\n")
        print(f"Current Input: {acfp}")
        print("Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:")
        user_input = input()
        # Check for exit conditions
        if user_input == "$":
            print(f"\n\nFinal Input:\n{acfp}\n\n")
            break
        elif user_input == "r" or user_input == "R":
            acfp = []
            print("Input cleared")
            continue
        
        if is_balanced_structure(user_input):

            # Check if the user input contains only ".", "(", and ")"
            
            if all(char == "." or char in ("(", ")") for char in user_input):
                if len(user_input) == len(acfp) + 1:
                    acfp.append(user_input)
                else:
                    print("Please add 1 character per step")
            else:
                print("Error: Invalid character in the folding path. Only '.', '(', and ')' are allowed.")
        else:
            print("Structure is not balanced -> closing/opening brackets don't match")

    return acfp

def couple(pair):
	if pair[0][0].upper() == pair[1][0].upper() and pair[0] != pair[1]:
		return True


def get_base_pairings_dict(dot_bracket):
    stack = []
    pairings_dict = {}

    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening_index = stack.pop()
                pairings_dict[opening_index] = i
                pairings_dict[i] = opening_index
    
    return pairings_dict


def acfp_to_domainfp(acfp,d_seq):
    """Takes an abstract folding path and a domain level sequence and converts it into a domain level path 
        Example: 
            acfp: [["."],["()"],[".()"],["()()"]],d_seq="b l b* a* l a b c d l d* c* b*"

            returns [[".."],["(.).."],["..((.))..."],["(.)...(((.)))"]]

        Args: 

            acfp(list): each sublist corresponds to one step example wrong not sublists
            d_seq(str): domain sequence
    """
    domain_fp = []
    d_seq = d_seq.split()


    #print(f'{d_seq = }')
    #print(f'{acfp =}')


    for x,cur_path in enumerate(acfp):
        module_index = 0

        base_pair_dict = get_base_pairings_dict(cur_path)
        path = ""
        stack = []
        #print(f'{cur_path = }')
        #print(' '.join(d_seq))
        struc_stack = []


        for i,domain in enumerate(d_seq):
            #print(f'{path = }, {domain = }, {cur_path[module_index]}')
            if domain[0] == "S":    
                module_index += 1

                if module_index == x + 1:
                    path += '.'
                    domain_fp.append(path)
                    break 

                else:
                    path += "."
                    

            elif cur_path[module_index] == ".":
                path += "."
                
            elif cur_path[module_index] == "(" and domain[0] != "S": 
                j = i
                k = module_index
                #print(f'{k = },{domain = }')
                added_notation = False
                while k < len(cur_path) and j <= len(d_seq) -1 :
                    if cur_path[k] == ")" and k == base_pair_dict[module_index]:
                        if couple((domain,d_seq[j])):
                            path += "("
                            struc_stack.append('(')
                            added_notation = True
                            stack.append(domain)
                            break

                    if cur_path[k] == "(" and k == base_pair_dict[module_index]:
                        if struc_stack:
                            if struc_stack[-1] == ")":
                                path += "."
                                added_notation = True

                            break

                    if d_seq[j][0] == "S":
                        k += 1 
                    j += 1

        
                if added_notation != True:
                    path += "."
                    

            elif cur_path[module_index] == ")" and domain[0] != "S":
                j = i 
                k = module_index
                added_notation = False 
                #print(f'{stack = }')
                if stack:
                    if couple((stack[-1],domain)) and struc_stack[-1] == "(":
                            path += ")"
                            stack.pop()
                            struc_stack.pop()
                    else:
                        path += "."

                else:
                        path += "."




    return domain_fp

def domainfp_to_acfp(domain_fp, d_seq):
    """Takes a domain level path and a domain level sequence and converts it into an abstract folding path.
        Example: 
            domain_fp: [[".."], ["(.).."], ["..((.))..."], ["(.)...(((.)))"]]
            d_seq: "b l b* a* l a b c d l d* c* b*"

            returns: [["."], ["()"], [".()"], ["()()"]]

        Args: 
            domain_fp(list): each sublist corresponds to one step
            d_seq(str): domain sequence
    """
    
    #print(f"{domain_fp = }")
    #print(f"{d_seq = }")

    acfp = []
    d_seq = d_seq.split()

    logic_indices = [index for index,entry in enumerate(d_seq) if entry.startswith("L") ]

    for i,step in enumerate(domain_fp):
        cur_path = []
        for struct in step:
            cur_struct = [struct[x] for x in logic_indices[0:i+1]]
            cur_struct = "".join(cur_struct)
            cur_path.append(cur_struct)
        
        cur_path = list(set(cur_path))
        acfp.append(cur_path)
    
    return acfp

def call_findpath(seq, ss1, ss2, md, fpw, mxb = float('inf')):
    """ Call ViennaRNA findpath. Modified from DrTransformerd_length
    """
    fc = RNA.fold_compound(seq)
    
    if mxb == float('inf'):
        path = fc.path_findpath(ss1, ss2, width = fpw)
        
    else:
        e1 = round(fc.eval_structure(ss1), 4)
        dcal_bound = int(round((mxb + e1)))
        path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = fpw)
        
    del fc

    if len(path):
        mypath = []
        barrier = None
        for step in path:
            struct = step.s
            energy = float(round(step.en,4))
            mypath.append((struct, energy))
            if barrier is None or barrier < energy:
                barrier = energy
        barrier -= float(round(path[0].en,4))
        del step, path # potentially a good idea.
        return mypath, barrier
    return None, None


if __name__=="__main__":
   

    print(' ')
    seq = "AAUACAUGCCCAAUCAUGUAUUUUAUCAUUUAUUGAAAUACAUGAUUCCUAAAAAACGCAAGGAAUCAUGUAUUUUAAUGCGCUCCUUCUCCCCUGCGUAUUGAAAUACAUGAUUCCUUGUAUCGCCACCCACCCCUAGAUGCAAGGAAUCAUGUAUUUUAAUGCGCAGGACCACUUUAUUUUUUC" 
    ss1 = '...........((((((((((((((........))))))))))))))..........((((((((((((((((((((((((((............)))))))))))))))))))))))))).................................................................' 
    ss2 = '((((((((......))))))))........((((((((((((((((((((..........))))))))))))))))))))............((((((((((((((((((((((((((((((((..............))))))))))))))))))))))))))))))))................'        
    mypath,barrier = call_findpath(seq,ss1,ss2,md=20,fpw = 20)
    for path in mypath:
        print(path)
    print(barrier)
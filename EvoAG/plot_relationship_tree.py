import graphviz
import os
import re
import math

def get_human_readable_size(size_bytes):
    """Converts a size in bytes to a human-readable format (e.g., KB, MB)."""
    if size_bytes == 0:
        return "0B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size_bytes, 1024)))
    p = math.pow(1024, i)
    s = round(size_bytes / p, 1)
    return f"{s}{size_name[i]}"

def parse_bact_id(file_path, files_dir=None):
    """
    Parses the bactID.txt file to extract relationships.
    Each line is expected in the format:
    child1+child2+...+childN OutsidE outgroup = parent
    """
    dot = graphviz.Digraph('BactRelationships', comment='Bacterial Relationship Graph')
    dot.attr('node', shape='box', style='rounded')
    dot.attr(rankdir='LR') # Left to Right layout

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            match = re.match(r'(.+?)\s+OutsidE\s+(.+?)\s+=\s+(.+)', line)
            if not match:
                print(f"Skipping malformed line: {line}")
                continue

            children_str, outgroup, parent = match.groups()
            
            parent = parent.strip()
            children = [child.strip() for child in children_str.split('+')]
            
            # --- Create nodes with updated labels ---
            all_node_names = {parent, *children, *outgroup.split('+')}
            for node_name_raw in all_node_names:
                node_name = node_name_raw.strip()
                if not node_name: continue

                label = node_name
                if files_dir:
                    # Check for the exact node name or with a .fasta extension
                    possible_filenames = [node_name, f"{node_name}.fasta"]
                    found_filepath = None
                    for filename in possible_filenames:
                        filepath = os.path.join(files_dir, filename)
                        if os.path.isfile(filepath):
                            found_filepath = filepath
                            break
                    
                    if found_filepath:
                        try:
                            size_bytes = os.path.getsize(found_filepath)
                            readable_size = get_human_readable_size(size_bytes)
                            label = f"{node_name}+{readable_size}"
                        except OSError as e:
                            print(f"Warning: Could not read size of file {found_filepath}: {e}")
                dot.node(node_name, label)
            
            # --- Create edges from parent to children ---
            for child in children:
                dot.edge(parent, child)

    return dot

def main():
    """
    Main function to generate and save the relationship graph.
    Configure input/output paths below instead of using CLI arguments.
    """
    input_file = "/home/wangyiwei/BactAGmuti/BactAG-output2026-02-11 05:26:24/bactID.txt"
    output_basename = "/home/wangyiwei/BactAGmuti/relationship_graph"
    files_dir = "/home/wangyiwei/BactAGmuti/BactAG-output2026-02-11 05:26:24/output"

    if not os.path.exists(input_file):
        print(f"Error: Input file not found at {input_file}")
        return

    if files_dir and not os.path.isdir(files_dir):
        print(f"Error: Specified files directory not found at {files_dir}")
        return

    print(f"Parsing file: {input_file}")
    if files_dir:
        print(f"Annotating with file sizes from: {files_dir}")

    graph = parse_bact_id(input_file, files_dir)

    output_path_png = f"{output_basename}.png"
    output_path_pdf = f"{output_basename}.pdf"

    try:
        print(f"Rendering graph to {output_path_png} and {output_path_pdf}...")
        graph.render(output_basename, format='png', view=False, cleanup=True)
        graph.render(output_basename, format='pdf', view=False, cleanup=True)
        print("Graph generation complete.")
        print(f"PNG saved to: {os.path.abspath(output_path_png)}")
        print(f"PDF saved to: {os.path.abspath(output_path_pdf)}")
    except Exception as e:
        print(f"An error occurred during graph rendering: {e}")
        print("Please ensure Graphviz is installed and in your system's PATH.")
        print("You can install it with: sudo apt-get install graphviz (on Debian/Ubuntu) or brew install graphviz (on macOS).")

if __name__ == '__main__':
    main()

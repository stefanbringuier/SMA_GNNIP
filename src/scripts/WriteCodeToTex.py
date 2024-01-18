import re
import sys

import paths


def extract_function(source_file, function_name):
    with open(source_file, "r") as file:
        code = file.read()

    pattern = re.compile(rf"(def {function_name}\(.*?:\n(?:\s+.*\n)+)", re.MULTILINE)
    match = pattern.search(code)

    if match:
        return match.group(1)
    else:
        raise ValueError(f"Function '{function_name}' not found in the file.")


def write_to_tex(output_file, function_code):
    """write extracted code in LaTeX listing environment"""
    tex_content = (
        f"\\begin{{lstlisting}}[language=Python]\n{function_code}\\end{{lstlisting}}\n"
    )

    with open(output_file, "w") as file:
        file.write(tex_content)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python script.py <source_python_file> <function_name> <output_tex_file>"
        )
        sys.exit(1)

    source_python_file = paths.scripts / sys.argv[1]
    function_name = sys.argv[2]
    output_tex_file = paths.output / sys.argv[3]

    try:
        python_function_code = extract_function(source_python_file, function_name)
        write_to_tex(output_tex_file, python_function_code)
        print(f"Function {function_name} has been written to {output_tex_file}")
    except ValueError as e:
        print(e)

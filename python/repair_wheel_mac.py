import argparse
import os
import shutil
import subprocess
import zipfile
from pathlib import Path
import glob


def extract_wheel(wheel_path, extract_to):
    """Extract the .whl file to a temporary directory."""
    with zipfile.ZipFile(wheel_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

def find_dylibs(so_path):
    """Find dynamic libraries linked by the shared object."""
    print('so_path', so_path)

    result = subprocess.run(['otool', '-L', so_path], stdout=subprocess.PIPE)
    dylibs = []
    for line in result.stdout.decode().splitlines():
        print(line)
        if '.dylib ' in line and not 'libSystem' in line:
            dylib = line.split(' ')[0].strip()
            if '@' in dylib:
                print('relative found. replace: ', dylib)
                dylib = dylib.replace('@rpath', os.path.dirname(so_path))
                print('result: ', dylib)

            dylibs.append(dylib)
    return dylibs

def copy_dylibs(dylibs, target_dir):
    """Copy dynamic libraries to the target directory."""
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    for dylib in dylibs:
        if dylib.startswith("@"):
            print(f"skip copying runtime-relative path: {dylib}")
            continue
        print('copy: ', dylib, target_dir)
        try:
            shutil.copy(dylib, target_dir)
        except PermissionError:
            pass

def fix_rpaths(so_path, lib_list):
    """Fix the rpaths in the shared object to point to the copied libraries."""
    for dylib_path in lib_list:
        dylib = os.path.basename(dylib_path)
        print('repair: ', 'install_name_tool', '-change', dylib_path, f'@loader_path/.dylibs/{dylib}', so_path)
        subprocess.run(['install_name_tool', '-change', dylib_path, f'@loader_path/.dylibs/{dylib}', so_path])
    # install_name_tool -change /usr/local/opt/gcc/lib/gcc/11/libgfortran.5.dylib @loader_path/../../../symgroupy.dylibs/libgfortran.5.dylib symgrouplib.cpython-311-darwin.so

def repair_wheel(wheel_path, output_dir):
    """Repair the wheel by bundling dynamic libraries."""

    extract_to = '/tmp/extracted_wheel'
    dylib_dir = os.path.join(extract_to, 'symgroupy/.dylibs')
    os.makedirs(dylib_dir, exist_ok=True)

    print('wheel_path: ', wheel_path, extract_to)
    extract_wheel(wheel_path, extract_to)

    fitxers = [f for f in os.listdir(extract_to)]
    print('fitxers: ', fitxers)

    #data_dir = glob.glob(os.path.join(extract_to, 'symgroupy-*.data'))[0]
    #so_lib_dir = os.path.join(data_dir, 'purelib', 'symgroupy')
    so_lib_dir = os.path.join(extract_to, 'symgroupy')

    for so_path in glob.glob(os.path.join(so_lib_dir, '*.so')):
        dylibs = find_dylibs(so_path)

        # second level of dependencies of dependencies
        second_level = []
        for lib in dylibs:
            second_level += find_dylibs(lib)

        dylibs += second_level

        copy_dylibs(dylibs, dylib_dir)
        fix_rpaths(so_path, dylibs)

    # Repack the repaired wheel
    #repaired_wheel_path = os.path.join(output_dir, os.path.basename(wheel_path))
    #with zipfile.ZipFile(repaired_wheel_path, 'w') as zip_ref:
    #    for root, dirs, files in os.walk(extract_to):
    #        for file in files:
    #            zip_ref.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), extract_to))

    # Repack the repaired wheel and regenerate RECORD
    subprocess.run(
        ["python", "-m", "wheel", "pack", extract_to, "--dest-dir", output_dir],
        check=True,
    )

def main():
    parser = argparse.ArgumentParser(description="Repair a Python wheel by bundling dynamic libraries.")
    parser.add_argument('wheel', type=str, help="Path to the .whl file.")
    parser.add_argument('-w', '--wheelhouse', type=str, default='wheelhouse', help="Directory to store the repaired wheel.")
    args = parser.parse_args()

    wheel_path = Path(args.wheel).resolve(strict=True)
    wheelhouse = Path(args.wheelhouse).resolve()
    wheelhouse.mkdir(parents=True, exist_ok=True)

    repair_wheel(wheel_path, wheelhouse)
    print(f"Repaired wheel saved to: {wheelhouse}")

if __name__ == "__main__":
    main()

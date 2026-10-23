"""Run anopheles_larva.py with the `bpy` Python module instead of the Blender app.
Usage: python run_with_bpy.py --out ./out [--print]
Requires: pip install bpy (each bpy release needs a specific Python version)
"""
import os, sys, runpy
sys.argv = ["blender", "--"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "anopheles_larva.py"), run_name="__main__")

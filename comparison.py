import os
from PIL import Image
import numpy as np

def load_ppm(path):
    return np.array(Image.open(path))

def normalized_distance(img1, img2):
    diff = np.abs(img1.astype(np.float32) - img2.astype(np.float32))
    return np.mean(diff) / 255.0

def main():
    main_folder = '.'
    render_folder = 'renders'
    
    # Match PPM files by name
    main_images = sorted(f for f in os.listdir(main_folder) if f.endswith('.ppm'))
    render_images = sorted(f for f in os.listdir(render_folder) if f.endswith('.ppm'))

    common = set(main_images).intersection(render_images)
    if not common:
        print("No matching PPM files found.")
        return

    distances = []
    for fname in sorted(common):
        img1 = load_ppm(os.path.join(main_folder, fname))
        img2 = load_ppm(os.path.join(render_folder, fname))
        if img1.shape != img2.shape:
            print(f"Skipping {fname}: image shapes do not match.")
            continue
        dist = normalized_distance(img1, img2)
        distances.append(dist)

    if not distances:
        print("No comparable images found.")
        return

    x = sum(distances) / len(distances)
    
    try:
        with open("render.cpp", "r") as f:
            code_length = len(f.read())
    except FileNotFoundError:
        print("render.cpp not found.")
        return

    final_score = code_length * (1 + x)
    
    print(f"Normalized distance (x): {x:.4f}")
    print(f"Code length (bytes): {code_length}")
    print(f"Final score: {final_score:.2f}")

if __name__ == "__main__":
    main()

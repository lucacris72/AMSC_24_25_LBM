import pandas as pd
import matplotlib.pyplot as plt
import imageio
import os

# Parametri
csv_files = ["test.csv"]  # List of CSV files
gif_filename = "simulation.gif"
image_folder = "images"

# Create folder to save data if doesn't exist
if not os.path.exists(image_folder):
    os.makedirs(image_folder)

# List to save name of generated images
images = []

# Read CSV files and create images
for i, csv_file in enumerate(csv_files):
    # Load data
    data = pd.read_csv(csv_file)

    # Create image with Matplotlib
    plt.figure(figsize=(8, 6))

    # Plot of density as heatmap
    plt.scatter(data['x'], data['y'], c=data['rho'], cmap='viridis', s=20)
    plt.colorbar(label='Densità (rho)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f"Passo temporale {i + 1}")

    # Salve image
    image_path = os.path.join(image_folder, f"frame_{i}.png")
    plt.savefig(image_path)
    images.append(image_path)

    # Close image to free memory
    plt.close()

# Create GIF from images
with imageio.get_writer(gif_filename, mode='I', duration=0.5) as writer:
    for image_path in images:
        image = imageio.imread(image_path)
        writer.append_data(image)

# Rimuove tmp images
for image_path in images:
    os.remove(image_path)

print(f"GIF salvata come {gif_filename}")

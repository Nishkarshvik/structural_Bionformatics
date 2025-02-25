#!/usr/bin/env python3

import os
import torch
import numpy as np
from PIL import Image
import torchvision.transforms as transforms
from realesrgan.utils import RealESRGANer
from basicsr.archs.rrdbnet_arch import RRDBNet  # Required for initializing the model

# Check for GPU availability
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Initialize the RRDBNet model
model = RRDBNet(num_in_ch=3, num_out_ch=3, num_feat=64, num_block=23, num_grow_ch=32, scale=4)

# Load Real-ESRGAN
upscaler = RealESRGANer(
    scale=4,
    model_path="https://github.com/xinntao/Real-ESRGAN/releases/download/v0.1.0/RealESRGAN_x4plus.pth",
    model=model,
    tile=0,  # Set tiling for large images, change if needed
    tile_pad=10,
    pre_pad=0,
    half=torch.cuda.is_available(),  # Use half-precision if CUDA is available
    device=device,
)

# Get list of all JPG files in the current directory
jpg_files = [f for f in os.listdir() if f.lower().endswith(".jpg")]

if not jpg_files:
    print("No JPG files found in the current directory.")
    exit()

# Process each JPG file
for jpg_file in jpg_files:
    try:
        # Open image
        img = Image.open(jpg_file).convert("RGB")

        # Upscale the image (Real-ESRGAN outputs a NumPy array, NOT a tensor)
        output_np, _ = upscaler.enhance(np.array(img), outscale=4)

        # Ensure pixel values are in the correct range
        output_np = np.clip(output_np, 0, 255).astype(np.uint8)

        # Convert to PIL Image
        upscaled_img = Image.fromarray(output_np)

        # Save as PNG with the same name
        output_file = f"{os.path.splitext(jpg_file)[0]}_upscaled.png"
        upscaled_img.save(output_file, format="PNG")

        print(f"Upscaled: {jpg_file} → {output_file}")

    except Exception as e:
        print(f"Error processing {jpg_file}: {e}")

print("Upscaling complete!")


import fitz  # PyMuPDF
import sys

def is_blank_area(page, rect, threshold=0.1):
    """Check if a rectangular area is blank by analyzing pixel data."""
    pix = page.get_pixmap(clip=rect, dpi=72)
    pixels = pix.samples
    # Consider area blank if most pixels are white (RGB: 255,255,255)
    white_pixels = sum(1 for i in range(0, len(pixels), 3) if pixels[i:i+3] == b'\xff\xff\xff')
    total_pixels = pix.width * pix.height
    return (white_pixels / total_pixels) > (1 - threshold)

def detect_content_bounds(page, margin=10):
    """Detect the bounding box of non-blank content on a page."""
    rect = page.rect  # Full page rectangle
    x0, y0, x1, y1 = rect.x0, rect.y0, rect.x1, rect.y1
    
    # Binary search to find content boundaries
    def find_edge(start, end, is_x, is_top_or_left):
        while abs(end - start) > 1:
            mid = (start + end) / 2
            if is_x:
                test_rect = fitz.Rect(mid, y0, mid + 1, y1) if is_top_or_left else fitz.Rect(x0, y0, mid, y1)
            else:
                test_rect = fitz.Rect(x0, mid, x1, mid + 1) if is_top_or_left else fitz.Rect(x0, y0, x1, mid)
            if is_blank_area(page, test_rect):
                if is_top_or_left:
                    start = mid
                else:
                    end = mid
            else:
                if is_top_or_left:
                    end = mid
                else:
                    start = mid
        return end if is_top_or_left else start

    # Find content boundaries
    left = find_edge(x0, x1, True, True)
    right = find_edge(x0, x1, True, False)
    top = find_edge(y0, y1, False, True)
    bottom = find_edge(y0, y1, False, False)

    # Add margin to avoid cutting content
    left = max(x0, left - margin)
    right = min(x1, right + margin)
    top = max(y0, top - margin)
    bottom = min(y1, bottom + margin)

    return fitz.Rect(left, top, right, bottom)

def crop_pdf(input_path, output_path, margin=10):
    """Crop blank space from all pages of a PDF and save to output_path."""
    try:
        doc = fitz.open(input_path)
        for page in doc:
            # Get content bounds
            content_rect = detect_content_bounds(page, margin)
            # Set new crop box
            page.set_cropbox(content_rect)
        
        # Save the cropped PDF
        doc.save(output_path, garbage=4, deflate=True)
        doc.close()
        print(f"Cropped PDF saved as {output_path}")
    except Exception as e:
        print(f"Error processing PDF: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python crop_pdf.py input.pdf output.pdf")
        sys.exit(1)
    input_pdf = sys.argv[1]
    output_pdf = sys.argv[2]
    crop_pdf(input_pdf, output_pdf)

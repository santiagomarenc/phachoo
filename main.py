import tkinter as tk
from tkinter import filedialog, messagebox

def on_button_click():
    file_path = filedialog.askopenfilename(
        title="Selecciona un archivo .tsv",
        filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
    )
    if file_path:
        if file_path.lower().endswith('.tsv'):
            messagebox.showinfo("Archivo seleccionado", f"Archivo válido: {file_path}")
            # Aquí puedes continuar con el análisis del archivo
        else:
            messagebox.showerror("Error", "Por favor selecciona un archivo con extensión .tsv")

def main():
    root = tk.Tk()
    root.title("Phachoo Data Analysis")
    root.geometry("1200x800")  # Set window size (optional)
    button = tk.Button(root, text="Seleccionar archivo .tsv", command=on_button_click)
    button.pack(pady=80)
    button = tk.Button

    root.mainloop()



if __name__ == "__main__":
    main()

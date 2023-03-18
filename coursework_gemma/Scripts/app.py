from flask import Flask, render_template
from dog_module import align_dog_breeds

app = Flask(__name__)

@app.route('/')
def index():
    test_sequence = '/home/gemma/PycharmProjects/pythonProject/coursework_gemma/Data_files/mystery.fa'
    best_breed, best_percentage_chance, p_value, doggy_dict = align_dog_breeds(test_sequence)
    return render_template('index.html', best_breed=best_breed, best_percentage_chance=best_percentage_chance, p_value=p_value, doggy_dict=doggy_dict)

if __name__ == "__main__":
    app.run(debug=True)

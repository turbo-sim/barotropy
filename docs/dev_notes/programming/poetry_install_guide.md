
## ✅ **Poetry installation steps (Git Bash)**

Use the following command block to install Poetry in Windows using Git Bash:

```bash
curl -sSL https://install.python-poetry.org | python -
echo '' >> ~/.bashrc
echo '# Add poetry to path' >> ~/.bashrc
echo 'export PATH="$PATH:$APPDATA\Python\Scripts"' >> ~/.bashrc
source ~/.bashrc
poetry --version
```


### 📝 **What each command does**

```bash
curl -sSL https://install.python-poetry.org | python -
```
- 📥 Downloads and runs the official Poetry installer using `curl`.
- 🐍 `python -` executes the downloaded script using Python.
- 🛠 Installs Poetry under `%APPDATA%\Python\Scripts`.

---

```bash
echo '' >> ~/.bashrc
echo '# Add poetry to path' >> ~/.bashrc
echo 'export PATH="$PATH:$APPDATA\Python\Scripts"' >> ~/.bashrc
```
- ✍️ Appends a line to your `.bashrc` to add Poetry's bin directory to your `PATH`.
- 🔁 Converts backslashes in `$APPDATA` to forward slashes (`/`), which Git Bash understands.
- ✅ Works regardless of the username or system.

---

```bash
source ~/.bashrc
```
- 🔁 Reloads your shell config so the updated `PATH` takes effect immediately.

---

```bash
poetry --version
```
- ✅ Verifies that Poetry is installed and correctly added to the path.

---



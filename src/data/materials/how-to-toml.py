import tomllib


if __name__ == "__main__":
    # Load the TOML file
    with open("src/data/materials/ss310.toml", "rb") as f:
        ss310 = tomllib.load(f)

    # Print the name and density from the loaded data
    print(ss310["name"])
    print(ss310["density"])

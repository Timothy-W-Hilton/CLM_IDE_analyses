from geopy.geocoders import Nominatim

if __name__ == "__main__":

    geolocator = Nominatim()
    location = geolocator.reverse("40.0, -100")
    print(location.address)

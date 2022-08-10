#include <iostream>
#include <thread>

using namespace std;

static bool finished = false;

void DoWork(int x)
{
    while (!finished)
    {
        cout << x << endl;
        this_thread::sleep_for(chrono::milliseconds(1000));
    }
}

int main() {

    thread worker(DoWork, 1);
    cin.get();
    finished = true;
    
    worker.join();
    cout << "Finished" << endl;
}
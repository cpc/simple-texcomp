package fi.tuni.cpc.texcompapp;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.EditText;
import android.widget.Spinner;
import android.widget.Switch;
import android.widget.TextView;

import java.io.File;

import fi.tuni.cpc.texcompapp.databinding.ActivityMainBinding;

public class MainActivity extends AppCompatActivity {

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

    private ActivityMainBinding binding;
    //private String inp_name = "/storage/35383430-2d36-4532-3000-2ba900000000/DCIM/OpenCamera/IMG_20201113_220310.jpg";
    private String inp_dir;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        binding = ActivityMainBinding.inflate(getLayoutInflater());
        setContentView(binding.getRoot());

        // Example of a call to a native method
        TextView tv = binding.sampleText;
        tv.setText("Hello World!");
    }

    public void setDir(View view) {
        EditText editText = (EditText) findViewById(R.id.editText);
        String dir_name = editText.getText().toString();

        File dir_path = new File(dir_name);

        //get the spinner from the xml.
        Spinner dropdown = findViewById(R.id.spinner1);
        //create a list of items for the spinner.
        //String[] items = new String[]{"1", "2", "three"};
        String[] items = dir_path.list();

        if (items == null) {
            TextView tv = binding.sampleText;
            tv.setText("No child paths found in a directory.");
            return;
        }

        //create an adapter to describe how the items are displayed, adapters are used in several places in android.
        //There are multiple variations of this, but this is the basic variant.
        ArrayAdapter<String> adapter = new ArrayAdapter<>(this, android.R.layout.simple_spinner_dropdown_item, items);
        //set the spinners adapter to the previously created one.
        dropdown.setAdapter(adapter);

        dropdown.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            @Override
            public void onItemSelected(AdapterView<?> parent, View view,
                                       int position, long id) {
                inp_dir = dir_name + "/" + (String) parent.getItemAtPosition(position);
                TextView tv = binding.sampleText2;
                tv.setText("dir: " + inp_dir);
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent) {
                // TODO Auto-generated method stub
            }
        });
    }

    public void runTranscoder(View view) {
        Switch jpeg_switch = (Switch) findViewById(R.id.switch1);
        String format;
        if (jpeg_switch.isChecked()) {
            format = "jpeg";
        } else {
            format = "astc";
        }

        //String resolution = getImgInfo(name);
        String exit_code = encodeDirectory(inp_dir, format);

        TextView tv = binding.sampleText;
        tv.setText("ret: " + exit_code);
    }

    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public native String getImgInfo(String inp_name);
    public native String encodeDirectory(String inp_name, String format);
}
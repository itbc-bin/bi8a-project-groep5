$(".download-button").on("click", () => {
    $.ajax({
        type: 'GET',
        url: '/download',
        contentType: 'csv',
        cache: false,
        processData: false,
        async: false,
    })
});

$(function () {
    $('#submit').click(() => {
        const form_data = new FormData($('#upload-form')[0]);
        $.ajax({
            type: 'POST',
            url: '/upload_file',
            data: form_data,
            contentType: false,
            processData: false,
            dataType: 'json'
        }).done((data) => {
            alert("file succesfully uploaded!")
        })
    });
});
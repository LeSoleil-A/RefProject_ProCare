// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2018 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#pragma once

#include <Core/Geometry/Image.h>
#include <Visualization/Shader/ShaderWrapper.h>

namespace open3d {

namespace glsl {

class ImageShader : public ShaderWrapper {
public:
    ~ImageShader() override { Release(); }

protected:
    ImageShader(const std::string &name) : ShaderWrapper(name) { Compile(); }

protected:
    bool Compile() final;
    void Release() final;
    bool BindGeometry(const Geometry &geometry,
                      const RenderOption &option,
                      const ViewControl &view) final;
    bool RenderGeometry(const Geometry &geometry,
                        const RenderOption &option,
                        const ViewControl &view) final;
    void UnbindGeometry() final;

protected:
    virtual bool PrepareRendering(const Geometry &geometry,
                                  const RenderOption &option,
                                  const ViewControl &view) = 0;
    virtual bool PrepareBinding(const Geometry &geometry,
                                const RenderOption &option,
                                const ViewControl &view,
                                Image &image) = 0;

protected:
    GLuint vertex_position_;
    GLuint vertex_position_buffer_;
    GLuint vertex_UV_;
    GLuint vertex_UV_buffer_;
    GLuint image_texture_;
    GLuint image_texture_buffer_;
    GLuint vertex_scale_;

    GLHelper::GLVector3f vertex_scale_data_;
};

class ImageShaderForImage : public ImageShader {
public:
    ImageShaderForImage() : ImageShader("ImageShaderForImage") {}

protected:
    virtual bool PrepareRendering(const Geometry &geometry,
                                  const RenderOption &option,
                                  const ViewControl &view) final;
    virtual bool PrepareBinding(const Geometry &geometry,
                                const RenderOption &option,
                                const ViewControl &view,
                                Image &render_image) final;
};

}  // namespace glsl

}  // namespace open3d
